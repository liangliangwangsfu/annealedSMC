package smcsampler;



import static nuts.util.CollUtils.map;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import nuts.math.HashGraph;
import nuts.math.Sampling;
import nuts.util.Arbre;
import pty.RootedTree;
import pty.RootedTree.RootingInfo;
import pty.RootedTree.Util.RootedTreeImpl;
import pty.UnrootedTree;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import goblin.Taxon;

public interface ProposalDistribution
{
	/**
	 * 
	 * @param current
	 * @param rand
	 * @return (proposedState, log proposal-ratio), where proposal-ratio =
	 *         Q(old|new)/Q(new|old)
	 */
	public Pair<UnrootedTree, Double> propose(UnrootedTree current, Random rand);

	public String description();

	public static class Options {
		@Option
		public double multiplicativeBranchProposalScaling = 2.0;
		@Option
		public boolean useMultiplicativeBranchProposal = true;
		@Option
		public boolean useGlobalMultiplicativeBranchProposal = true;
		@Option
		public boolean useStochasticNearestNeighborInterchangeProposal = true;
		@Option
		public boolean useStochasticNearestNeighborInterchangeProposalWithNbrsResampling = true;
		@Option
		public boolean useSubtreePruningRegraftingProposal = true;
		@Option
		public boolean useIndepBranchProp = false;

	}

	public static class Util {
		public static final ProposalDistribution.Options _defaultProposalDistributionOptions = new ProposalDistribution.Options();

		public static List<ProposalDistribution> proposalList(Options options,
				UnrootedTree nct, Random rand) {
			List<ProposalDistribution> result = new ArrayList<ProposalDistribution>();
			if (options.useMultiplicativeBranchProposal)
				result.add(new MultiplicativeBranchProposal(
						options.multiplicativeBranchProposalScaling, false));
			if (options.useGlobalMultiplicativeBranchProposal)
				result.add(new MultiplicativeBranchProposal(
						options.multiplicativeBranchProposalScaling, true));
			if (options.useStochasticNearestNeighborInterchangeProposal)
				result.add(new StochasticNearestNeighborInterchangeProposal());
			if (options.useStochasticNearestNeighborInterchangeProposalWithNbrsResampling)
				result.add(new StochasticNearestNeighborInterchangeProposal(
						true, options.multiplicativeBranchProposalScaling));
			if (options.useSubtreePruningRegraftingProposal)
				result.add(new SubtreePruningRegraftingProposal());
			if (options.useIndepBranchProp)
				result.add(new IndepBranchProposal());
			Collections.shuffle(result, rand);
			return result;
		}
	}

	public static class IndepBranchProposal implements ProposalDistribution {

		@Override
		public String description() {
			return "indepBL";
		}

		@Override
		public Pair<UnrootedTree, Double> propose(UnrootedTree current,
				Random rand) {
			UnorderedPair<Taxon, Taxon> edge = current.randomEdge(rand);
			double oldVal = current.branchLength(edge);
			double newVal = Sampling.sampleExponential(rand, 1.0);
			double logRatio = Sampling.exponentialLogDensity(1.0, oldVal)
					- Sampling.exponentialLogDensity(1.0, newVal);
			final UnrootedTree proposedTree = current.branchLengthNeighbor(
					edge, newVal);
			return Pair.makePair(proposedTree, logRatio);
		}

	}

	public static class MultiplicativeBranchProposal implements
			ProposalDistribution {
		private final double a; // higher will have lower accept rate
		private final boolean global;
		public UnorderedPair<Taxon, Taxon> selectedEdge = null;

		public MultiplicativeBranchProposal(double a, boolean global) {
			if (a <= 1)
				throw new RuntimeException();
			this.global = global;
			this.a = a;
		}

		public Pair<UnrootedTree, Double> propose(UnrootedTree current,
				Random rand) {
			double lambda=2*Math.log(a);
			double rvUnif = Sampling.nextDouble(rand, 0, 1);
			double m  = Math.exp(lambda*(rvUnif-0.5));			
			if (global) {
				double sum = 0.0;
				UnrootedTree proposedTree = current;				
				for (UnorderedPair<Taxon, Taxon> edge : current.edges()) {
					rvUnif = Sampling.nextDouble(rand, 0, 1);
					m  = Math.exp(lambda*(rvUnif-0.5));								
					Pair<UnrootedTree, Double> p = propose(proposedTree, rand,
							edge, m);					
					proposedTree = p.getFirst();
					sum += p.getSecond();					
				}
				return Pair.makePair(proposedTree, sum);
			} else if (selectedEdge != null){
				return propose(current, rand, selectedEdge, m);
			}
			return propose(current, rand, current.randomEdge(rand), m);
		}

		public static Pair<UnrootedTree, Double> propose(UnrootedTree current,
				Random rand, UnorderedPair<Taxon, Taxon> edge, double m) {
			final double newBL = m * current.branchLength(edge);
			final UnrootedTree proposedTree = current.branchLengthNeighbor(
					edge, newBL);
			// final NonClockTreeState proposedState =
			// current.copyAndChange(proposedTree);
             return Pair.makePair(proposedTree, Math.log(m));						
		}

		public String description() {
			return (global ? "g" : "") + "MB(" + a + ")";
		}
	}

	public static class StochasticNearestNeighborInterchangeProposal implements
			ProposalDistribution {
		private final boolean resampleNbrEdges;
		private final double a;
		public UnorderedPair<Taxon, Taxon> selectedEdge = null;

		public StochasticNearestNeighborInterchangeProposal() {
			this.resampleNbrEdges = false;
			this.a = -1;
		}

		public StochasticNearestNeighborInterchangeProposal(
				boolean resampleNbrEdges, double a) {
			this.resampleNbrEdges = resampleNbrEdges;
			this.a = a;
		}

		private UnorderedPair<Taxon, Taxon> lastEdge = null;

		public UnorderedPair<Taxon, Taxon> getLastEdgePicked() {
			return lastEdge;
		}

		public Pair<UnrootedTree, Double> propose(UnrootedTree current,
				Random rand) {
			// final UnorderedPair<Taxon,Taxon> edge =
			// current.randomNonTerminalEdge(rand);
			UnorderedPair<Taxon, Taxon> edge = null;
			if (selectedEdge != null)
				edge = selectedEdge;
			else
				edge = current.randomNonTerminalEdge(rand);
			// LogInfo.logsForce("Picked edge:" + edge);
			if (edge == null)
				return null;
			lastEdge = edge;
			// // one third probability to stay here (to make sure
			// aperiodic---although not that important if we are summing)
			// if (rand.nextDouble() < 1.0/3.0) return Pair.makePair(current,
			// 1.0);
			// else, pick one of the 2 neighbors at random
			final int picked = rand.nextInt(2);
			// LogInfo.logsForce("Picked nhbr:" + picked);
			UnrootedTree proposedTree = current.topologicalNeighbors(edge).get(
					picked); // NOTE: this could be made faster
			// NonClockTreeState proposedState =
			// current.copyAndChange(proposedTree);
			double sum = 0.0;
			if (resampleNbrEdges)
				for (UnorderedPair<Taxon, Taxon> nbrEdge : proposedTree
						.nbrEdges(edge)) {
					//final double m = Sampling.nextDouble(rand, 1.0 / a, a);
					double lambda=2*Math.log(a);
					double rvUnif = Sampling.nextDouble(rand, 0, 1);
					double m  = Math.exp(lambda*(rvUnif-0.5));	
					Pair<UnrootedTree, Double> p = MultiplicativeBranchProposal
							.propose(proposedTree, rand, nbrEdge, m);
					proposedTree = p.getFirst();
					sum += p.getSecond();
				}
			return Pair.makePair(proposedTree, sum);
		}

		public String description() {
			return "sNNI" + (resampleNbrEdges ? "+MB(" + a + ")" : "");
		}
	}

	public static class SubtreePruningRegraftingProposal implements
			ProposalDistribution {
		private final boolean resampleNbrEdges;
		private final double a;
		public UnorderedPair<Taxon, Taxon> selectedEdge = null;

		public SubtreePruningRegraftingProposal() {
			this.resampleNbrEdges = false;
			this.a = -1;
		}

		public SubtreePruningRegraftingProposal(
				boolean resampleNbrEdges, double a) {
			this.resampleNbrEdges = resampleNbrEdges;
			this.a = a;
		}

		private UnorderedPair<Taxon, Taxon> lastEdge = null;

		public UnorderedPair<Taxon, Taxon> getLastEdgePicked() {
			return lastEdge;
		}

		public Pair<UnrootedTree, Double> propose(UnrootedTree current,
				Random rand) {
			// System.out.println(current.toString());
			UnorderedPair<Taxon, Taxon> edge = null;
			if (selectedEdge != null)
				edge = selectedEdge;
			else
				edge = current.randomNonTerminalEdge(rand);
			// LogInfo.logsForce("Picked edge:" + edge);
			if (edge == null)
				return null;
			lastEdge = edge;
			// System.out.println(edge.getFirst().toString() + " "
			// + edge.getSecond().toString());

			double edgeBl=current.branchLength(edge);
			
			Pair<RootedTree, RootedTree> twoRootedTrees = divideOneUnrootedTree2twoRootedTrees(
					current, edge);

			// System.out.println("First rooted tree: "
			// + twoRootedTrees.getFirst());

			RootedTree secondTree = twoRootedTrees.getSecond();
			// System.out.println("Second rooted tree: " +
			// secondTree.toString());

			UnrootedTree firstTree = UnrootedTree.fromRooted(twoRootedTrees
					.getFirst());
			// System.out.println("First unrooted tree: " +
			// firstTree.toString());

			List<UnorderedPair<Taxon, Taxon>> allEdges=firstTree.edges();
			// allEdges.remove(edge);
			final int picked = rand.nextInt(allEdges.size());
			// System.out.println(allEdges.size() + " " + picked);
			
			UnorderedPair<Taxon, Taxon> regraftEdge = allEdges.get(picked);
			// System.out.println(regraftEdge);

			RootingInfo rooting = new RootingInfo(regraftEdge.getFirst(),
					regraftEdge.getSecond(), edge.getFirst(), rand.nextDouble());
			RootedTree rt = firstTree.reRoot(rooting);
			
			// System.out.println(rt.toString());
			Taxon regraftNewRoot =new Taxon("regraftNewRoot");
			
			@SuppressWarnings("unchecked")
			Arbre<Taxon> newArbre = Arbre.arbreWithChildren(
					regraftNewRoot,
					rt.topology(), secondTree.topology());
			 List<Arbre<Taxon>> newtreeRootChildren=newArbre.root().getChildren();			
			Map<Taxon, Double> blMap = map();
			blMap.putAll(rt.branchLengths());
			blMap.putAll(secondTree.branchLengths());
			blMap.put(newtreeRootChildren.get(0).getContents(), edgeBl*0.5);
			blMap.put(newtreeRootChildren.get(1).getContents(), edgeBl*0.5);
			RootedTree regraftRootedTree = new RootedTreeImpl(newArbre, blMap);
			UnrootedTree proposedTree = UnrootedTree
					.fromRooted(regraftRootedTree);
			double sum = 0.0;
			return Pair.makePair(proposedTree, sum);
		}

		public static UnrootedTree regraft(RootedTree c) {
			// Graph<Language> topo = new
			// Graph.HashGraph<Language>(Arbre.arbre2Tree(c.topology()));
			Map<UnorderedPair<Taxon, Taxon>, Double> branchLengths = new HashMap<UnorderedPair<Taxon, Taxon>, Double>();
			Set<Taxon> languages = new HashSet<Taxon>(c.topology()
					.nodeContents());

			Set<UnorderedPair<Taxon, Taxon>> edges = new HashSet<UnorderedPair<Taxon, Taxon>>();
			for (Arbre<Taxon> node : c.topology().nodes())
				if (!node.isRoot()) {
					UnorderedPair<Taxon, Taxon> edge = null;
					double bl = -1;
					if (node.getParent().isRoot()
							&& node.getParent().getChildren().size() == 2) { // special
																				// case
																				// if
																				// we
																				// are
																				// at
																				// the
																				// root:
																				// connect
																				// that
																				// two
																				// children,
																				// marginalizing
																				// root
																				// node
						languages.remove(c.topology().getContents());
						List<Arbre<Taxon>> parentChildren = node.getParent()
								.getChildren();
						edge = new UnorderedPair<Taxon, Taxon>(parentChildren
								.get(0).getContents(), parentChildren.get(1)
								.getContents());
						bl = c.branchLengths().get(
								parentChildren.get(0).getContents())
								+ c.branchLengths().get(
										parentChildren.get(1).getContents());
					} else {
						edge = new UnorderedPair<Taxon, Taxon>(
								node.getContents(), node.getParent()
										.getContents());
						bl = c.branchLengths().get(node.getContents());
					}
					edges.add(edge);
					branchLengths.put(edge, bl);
				}
			return new UnrootedTree(new HashGraph<Taxon>(languages, edges),
					branchLengths);
		}

		// public Pair<UnrootedTree, UnrootedTree>
		// divideOneUnrootedTree2twoUnrootedTrees(
		// UnrootedTree urt, UnorderedPair<Taxon, Taxon> selectEdge) {
		// Map<UnorderedPair<Taxon, Taxon>, Double> branchLengths = new
		// HashMap<UnorderedPair<Taxon, Taxon>, Double>();
		// Set<Taxon> vertices = new HashSet<Taxon>(urt.getTopology()
		// .vertexSet());
		// Set<UnorderedPair<Taxon, Taxon>> edges = new
		// HashSet<UnorderedPair<Taxon, Taxon>>();
		// new UnrootedTree(new HashGraph<Taxon>(vertices, edges),
		// branchLengths);
		//
		// urt.inducedBiPartitions2BranchMap();
		//
		// List<Taxon> allLeaves = urt.leaves();
		// Arbre<Taxon> topology =
		// Arbre.tree2Arbre(Graphs.toTree(urt.getTopology(),
		// selectEdge.getFirst()));
		//
		// topology.
		//
		//
		//
		//
		// Map<Arbre<Taxon>, Set<Taxon>> leavesMap = Arbre.leavesMap(topology);
		// for (Arbre<Taxon> key : leavesMap.keySet())
		// if (!key.isRoot())
		// {
		// Set<Taxon> clade = leavesMap.get(key);
		// double bl = branchLength(key.getContents(),
		// key.getParent().getContents());
		// Set<Taxon> complement = CollUtils.set(allLeaves);
		// complement.removeAll(clade);
		// result.setCount(new UnorderedPair<Set<Taxon>,Set<Taxon>>(complement,
		// clade), bl);
		// }
		// return result;
		//
		//
		//
		//
		//
		// // UnorderedPair<Taxon,Taxon> selectEdge = urt.randomEdge(rand);
		// RootingInfo rooting = new RootingInfo(selectEdge.getFirst(),
		// selectEdge.getSecond(), Taxon.dummy, 0.5);
		// RootedTree rt = urt.reRoot(rooting);
		// Arbre<Taxon> topo = rt.topology();
		// List<Arbre<Taxon>> rootChildren = topo.root().getChildren();
		// Arbre<Taxon> trtopo0 = rootChildren.get(0).copy(), trtopo1 =
		// rootChildren
		// .get(1).copy();
		//
		//
		// // int sz0 = trtopo0.getLeaves().size(), sz1 = trtopo1.getLeaves()
		// // .size();
		// Map<Taxon, Double> allbl = rt.branchLengths();
		// Map<Taxon, Double> subbl0 = map(), subbl1 = map();
		// for (Taxon tax : trtopo0.nodeContents())
		// subbl0.put(tax, allbl.get(tax));
		// for (Taxon tax : trtopo1.nodeContents())
		// subbl1.put(tax, allbl.get(tax));
		//
		// RootedTree subtr0 = new RootedTreeImpl(trtopo0, subbl0), subtr1 = new
		// RootedTreeImpl(
		// trtopo1, subbl1);
		// // int sz0 = subbl0.size(), sz1 = subbl1.size();
		// // if (sz0 <= sz1)
		// System.out.println(subtr0.toString());
		// System.out.println(subtr1.toString());
		// return Pair.makePair(subtr0, subtr1);
		// // else
		// // return Pair.makePair(subtr1, subtr0);
		// }

		public Pair<RootedTree, RootedTree> divideOneUnrootedTree2twoRootedTrees(
				UnrootedTree urt, UnorderedPair<Taxon, Taxon> selectEdge) {
			// UnorderedPair<Taxon,Taxon> selectEdge = urt.randomEdge(rand);
			RootingInfo rooting = new RootingInfo(selectEdge.getFirst(),
					selectEdge.getSecond(), Taxon.dummy, 0.5);
			RootedTree rt = urt.reRoot(rooting);
			Arbre<Taxon> topo = rt.topology();
			List<Arbre<Taxon>> rootChildren = topo.root().getChildren();

			Arbre<Taxon> trtopo0 = rootChildren.get(0).copy(), trtopo1 = rootChildren
					.get(1).copy(); // TODO: there is problem with this.

			// System.out.println(trtopo0.deepToString());
			// System.out.println(trtopo0.deepToString());
			// System.out.println(trtopo0.isRoot());
			// System.out.println(trtopo1.isRoot());
			// int sz0 = trtopo0.getLeaves().size(), sz1 = trtopo1.getLeaves()
			// .size();
			Map<Taxon, Double> allbl = rt.branchLengths();
			Map<Taxon, Double> subbl0 = map(), subbl1 = map();

			for (Taxon tax : trtopo0.nodeContents())
				if (!tax.equals(trtopo0.getContents()))
					subbl0.put(tax, allbl.get(tax));
			for (Taxon tax : trtopo1.nodeContents())
				if (!tax.equals(trtopo1.getContents()))
					subbl1.put(tax, allbl.get(tax));

			RootedTree subtr0 = new RootedTreeImpl(trtopo0, subbl0), subtr1 = new RootedTreeImpl(
					trtopo1, subbl1);

			// System.out.println(subbl0.keySet());

			// int sz0 = subbl0.size(), sz1 = subbl1.size();
			// if (sz0 <= sz1)
			// System.out.println(subtr0.toString());
			// System.out.println(subtr1.toString());
			return Pair.makePair(subtr0, subtr1);
			// else
			// return Pair.makePair(subtr1, subtr0);
		}


		public String description() {
			return "sSPR" + (resampleNbrEdges ? "+MB(" + a + ")" : "");
		}
	}

}
