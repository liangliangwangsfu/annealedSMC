package smcsampler;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Lists;

import conifer.trees.StandardNonClockPriorDensity;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.mcmc.PhyloSampler;
import pty.mcmc.PhyloSampler.NonClockTreePrior;
import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import ev.ex.TreeGenerators;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import fig.prob.Distrib;
import fig.prob.Gamma;
import goblin.Taxon;
import nuts.math.Graph;
import nuts.math.HashGraph;

public class AnnealingKernel implements SMCSamplerKernel<UnrootedTreeState>
{
	@Option public static double AnnealDeltaProposalRate = 10.0;    
	@Option
	public static int nAnnealing = 500;
	@Option public static boolean printBranchLengthMagnitudes = false;
	private final UnrootedTreeState initial;
	private double temperature = 0;
	private double newtemperature = 0;
	private double defaultTemperatureDifference =0;
	private boolean initializing = true;
	private int currentIter=0;

	private LinkedList<ProposalDistribution> proposalDistributions = null;
	private ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;

	public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	private double temperatureDifference = 0;

	public double getTemperatureDifference() {
		return temperatureDifference;
	}
	public void setTemperatureDifference(double temperatureDifference) {
		this.temperatureDifference = temperatureDifference;
		newtemperature = Math.min(temperature + temperatureDifference,
				1.0);
		//	LogInfo.logs("temperature: " + temperature +"  newtemperature: " + newtemperature);
		//	System.out.println("temperature: " + temperature +"  newtemperature: " + newtemperature);
	}
	public AnnealingKernel(UnrootedTreeState initial, double defaultTemperatureDifference, LinkedList<ProposalDistribution> proposalDistributions, ProposalDistribution.Options proposalOptions) 
	{
		this.initial = initial;		 
		this.proposalDistributions=proposalDistributions;
		this.proposalOptions=proposalOptions;
		this.defaultTemperatureDifference=defaultTemperatureDifference;
	}
	public UnrootedTreeState getInitial() { return initial; }

	public boolean isLastIter()
	{
		return temperature == 1.0;
	}

	public UnrootedTreeState sampleFromPrior(Random rand, UnrootedTreeState current)
	{
		//NonClockTreePrior; 
		Gamma exponentialPrior = Gamma.exponential(AnnealDeltaProposalRate);
//		RootedTree proprosedRTree = TreeGenerators.sampleExpNonclock(rand, current.getUnrootedTree().leaves(), AnnealDeltaProposalRate);
//		UnrootedTreeState proposedState = current.copyAndChange(UnrootedTree.fromRooted(proprosedRTree));
//		RootedTree proprosedRTree = generate(rand, exponentialPrior,current.getUnrootedTree().leaves()); 		
		UnrootedTreeState proposedState = current.copyAndChange(generate(rand, exponentialPrior,current.getUnrootedTree().leaves()));				
		return proposedState;
	}

	@Override
	public Pair<UnrootedTreeState, Double> next(
			Random rand,
			UnrootedTreeState current)
	{
		if (currentIter==0)  current = sampleFromPrior(rand, current);					
		//		PhyloSampler  sampler = new PhyloSampler();
		//		sampler.init(current);
		ProposalDistribution proposal = null; 
		UnrootedTreeState proposedState=null; 		 		
		while(proposal==null)
		{
			proposal = nextProposal(rand);
			Pair<UnrootedTree,Double> result = proposal.propose(current.getNonClockTree(), rand);
			if (result != null) // might happen e.g. when trying to do nni with 3 leaves
			{
				proposedState = current.copyAndChange(result.getFirst());
				final double logProposalRatio = result.getSecond();
				double currentPrior = current.getLogPrior(), proposePrior = proposedState
						.getLogPrior();
				double logLikRatio = proposePrior+newtemperature* proposedState.getLogLikelihood()
				- (currentPrior+newtemperature* current.getLogLikelihood());
				final double ratio = Math.min(1,
						Math.exp(logProposalRatio + logLikRatio));
				if (Double.isNaN(ratio))
					throw new RuntimeException();
				if (rand.nextDouble() >= ratio) {
					proposedState = current;
				}
			}
		}
		temperature = newtemperature;
		double logw = temperatureDifference * current.logLikelihood();
		return Pair.makePair(proposedState, logw);    
	}


	//	@Override
	//	public Pair<UnrootedTreeState, Double> next(Random rand,
	//			UnrootedTreeState current)
	//	{
	//		return (Pair) _next(rand, current, false);
	//	}

	//	@Override
	//	public int nIterationsLeft(UnrootedTreeState partialState) {	
	//		return (nAnnealing);
	//	}

	public void setCurrentIter(int  currentIter) {	
		this.currentIter=currentIter; 
	}

	public int getCurrentIter() {
		return this.currentIter;
	}


	private ProposalDistribution nextProposal(Random rand) {
		if (proposalDistributions.isEmpty())
			proposalDistributions.addAll(ProposalDistribution.Util
					.proposalList(proposalOptions, initial.getNonClockTree(),
							rand));
		return proposalDistributions.get(rand.nextInt(proposalDistributions
				.size()));
	}

	//	private double prior(UnrootedTree urt) {
	//		double result = 0.0;
	//		for (UnorderedPair<Taxon, Taxon> edge : urt.edges()) {
	//			result += Sampling.exponentialLogDensity(
	//					1.0 / AnnealDeltaProposalRate, urt.branchLength(edge));
	//		}
	//		return result;
	//	}

	public boolean isInitializing() {
		return initializing;
	}

	public void setInitializing(boolean initializing) {
		this.initializing = initializing;
	}

	public double getDefaultTemperatureDifference() {
		return defaultTemperatureDifference;
	}

	public void setDefaultTemperatureDifference(double defaultTemperatureDifference) {
		this.defaultTemperatureDifference = defaultTemperatureDifference;
	}

	
	/**
	   * Generate a tree with a topology uniformly distributed
	   * over bifurcating unrooted tree, using:
	   * 
	   *   The generation of random, binary unordered trees
	   *   George W. Furnas
	   *     see 2.1.2 (p.204)
	   */
	
	  public static UnrootedTree generate(
		      Random random, 
		      Distrib <Double> branchDistribution,
		      List<Taxon> leaves)
		  {
		    List<UnorderedPair<Taxon, Taxon>> result = new ArrayList<UnorderedPair<Taxon, Taxon>>();
		    List<Taxon> shuffled = Lists.newArrayList(leaves);
		    Collections.shuffle(shuffled, random);
		    Queue<Taxon> queue = Lists.newLinkedList(shuffled);
		    if (queue.isEmpty())
		      new RuntimeException();
		    
		    Taxon leaf1 = queue.poll();
		    if (queue.isEmpty())
		    	new RuntimeException();
		    
		    Taxon leaf2 = queue.poll();
		    result.add(new UnorderedPair<Taxon, Taxon>(leaf1, leaf2)); 
		    
		    int i=1;
		    while (!queue.isEmpty())
		    {
		      // pick a random edge
//		    	SamplingUtils.uniformFromCollection(random, result.getTopology().edgeSet());
		      UnorderedPair<Taxon, Taxon> edge =result.get(random.nextInt(result.size()));
		      Taxon internal=new Taxon("internal_" + (i++));		      
//		      TreeNode internal = TreeNode.nextUnlabelled();
		      Taxon newLeaf = queue.poll();
		      result.remove(edge);
		      result.add(new UnorderedPair<Taxon, Taxon>(newLeaf, internal));
		      result.add(new UnorderedPair<Taxon, Taxon>(internal, edge.getFirst()));
		      result.add(new UnorderedPair<Taxon, Taxon>(internal, edge.getSecond()));
		    }
		    		    
		    Graph<Taxon> topo=new HashGraph(new HashSet<UnorderedPair<Taxon, Taxon>>(result));
			Map<UnorderedPair<Taxon, Taxon>, Double> branchLengths = new HashMap<UnorderedPair<Taxon,Taxon>,Double>();
			
		    for (UnorderedPair<Taxon, Taxon> edge : result)		    			    
		    	branchLengths.put(edge, branchDistribution.sampleObject(random));
		    return new UnrootedTree(topo, branchLengths);
		  }
	
}