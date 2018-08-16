package evolmodel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import com.google.common.collect.Lists;
import conifer.trees.StandardNonClockPriorDensity;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import pty.smc.models.CTMC;
import smcsampler.SMCSamplerKernel;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import fig.prob.Distrib;
import fig.prob.Gamma;
import goblin.Taxon;
import nuts.math.Graph;
import nuts.math.HashGraph;

public class AnnealingKernelTreeEvolPara implements SMCSamplerKernel<UnrootedTreeEvolParameterState>
{
	@Option public static double AnnealDeltaProposalRate = 10.0;    
	@Option public static boolean printBranchLengthMagnitudes = false;
	private final UnrootedTreeEvolParameterState initial;
	private double temperature = 0;
	private double newtemperature = 0;
	private double defaultTemperatureDifference = 0;
	private boolean initializing = true;
	private int currentIter=0;
	private EvolutionModel evolModel = EvolutionModel.K2P; 
	private final Dataset dataset;
	private final StandardNonClockPriorDensity treePriorDensity; 

	private LinkedList<ProposalDistribution> proposalDistributions = null;
	private ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;

	private LinkedList<EvolutionParameterProposalDistribution> evolProposalDistributions = null;
	private EvolutionParameterProposalDistribution.Options evolProposalOptions = EvolutionParameterProposalDistribution.Util._defaultProposalDistributionOptions;

	
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
		double newTemp=temperature + temperatureDifference;
		if(newTemp>=1.0) 
		{
			newtemperature=1.0;
			this.temperatureDifference=1.0-temperature;			
		}
		else
		{
			newtemperature=newTemp;
			this.temperatureDifference=temperatureDifference;
		}			
		//		LogInfo.logs("temperature: " + temperature +"  newtemperature: " + newtemperature);
	}

	public AnnealingKernelTreeEvolPara(Dataset dataset, StandardNonClockPriorDensity treePriorDensity, UnrootedTreeEvolParameterState initial, double defaultTemperatureDifference, LinkedList<ProposalDistribution> proposalDistributions, ProposalDistribution.Options proposalOptions, LinkedList<EvolutionParameterProposalDistribution> evolProposalDistributions, EvolutionParameterProposalDistribution.Options evolProposalOptions) 
	{
		this.dataset = dataset;
		this.treePriorDensity = treePriorDensity;  
		this.initial = initial;		 
		this.proposalDistributions=proposalDistributions;
		this.proposalOptions=proposalOptions;
		this.defaultTemperatureDifference=defaultTemperatureDifference;
		this.evolProposalDistributions=evolProposalDistributions;
		this.evolProposalOptions=evolProposalOptions;		
	}
	
	public UnrootedTreeEvolParameterState getInitial() { return initial; }

	public boolean isLastIter()
	{
		return temperature >= 1.0;
	}

	public UnrootedTreeEvolParameterState sampleFromPrior(Random rand, UnrootedTreeEvolParameterState current)
	{ 
		Gamma exponentialPrior = Gamma.exponential(AnnealDeltaProposalRate);
		UnrootedTreeState proposedState = current.getUnrootedTreeState().copyAndChange(generate(rand, exponentialPrior, current.getUnrootedTreeState().getUnrootedTree().leaves()));
		current.getEvolParameter().sampleFromPrior(rand);		
		return new UnrootedTreeEvolParameterState(proposedState, current.getEvolParameter());
	}

	@Override
	public Pair<UnrootedTreeEvolParameterState, Double> next(
			Random rand,
			UnrootedTreeEvolParameterState current)
	{
		if (currentIter==0) 
		{
			current = sampleFromPrior(rand, current);
			return Pair.makePair(current, 0.0); 
		}
		ProposalDistribution proposal = nextProposal(rand);         
		UnrootedTree currenturt=current.getUnrootedTreeState().getNonClockTree();
		Pair<UnrootedTree,Double> result = proposal.propose(currenturt, rand);		
		UnrootedTreeState  proposedUrootedTreeState = null;				
		UnrootedTreeEvolParameterState  proposedState = null;
		double logw = temperatureDifference * current.getUnrootedTreeState().logLikelihood();
		if (result != null) // might happen e.g. when trying to do nni with 3 leaves
		{
			double logTargetDenCurrent=logTargetDensity(newtemperature,current.getUnrootedTreeState());
			EvolutionParameterProposalDistribution evolParaProposal = nextEvolParaProposal(rand);  
			Pair<EvolutionParameters, Double> evolProposalRe = evolParaProposal.propose(current.getEvolParameter(), rand);
			EvolutionParameters proposedEvolPar = evolProposalRe.getFirst();    // proposed evolutionary parameters
			CTMC ctmc = evolModel.instantiateCTMC(proposedEvolPar, dataset.nSites());			
			proposedUrootedTreeState = UnrootedTreeState.initFastState(result.getFirst(), dataset, ctmc, treePriorDensity);  // result.getFirst() gives the proposed tree and ctmc is updated with the newly proposed evolutionary parameters			
			final double logProposalRatio = result.getSecond()+evolProposalRe.getSecond();  // proposal ratio (tree and evolutionary parameters)
			double logLikRatio = logTargetDensity(newtemperature,proposedUrootedTreeState) -logTargetDenCurrent;  	//TODO: double check if logTargetDensity needs to be updated or not: i.e. adding the priors for the evolutionary parameters? 			
			final double ratio = Math.min(1, Math.exp(logProposalRatio + logLikRatio));
			if (Double.isNaN(ratio))
				throw new RuntimeException();
			if (rand.nextDouble() >= ratio) {
				proposedState = current;
			}else
			{
				proposedState = new UnrootedTreeEvolParameterState(proposedUrootedTreeState, proposedEvolPar);
			}
		}
		temperature = newtemperature;
		return Pair.makePair(proposedState, logw);    
	}

	public double logTargetDensity(double temperature, UnrootedTreeState uts)
	{	
		return uts.getLogPrior()+ temperature*uts.getLogLikelihood();
	}


	public void setCurrentIter(int  currentIter) {	
		this.currentIter=currentIter; 
	}

	public int getCurrentIter() {
		return this.currentIter;
	}


	private ProposalDistribution nextProposal(Random rand) {
		if (proposalDistributions.isEmpty())
			proposalDistributions.addAll(ProposalDistribution.Util
					.proposalList(proposalOptions, initial.getUnrootedTreeState().getNonClockTree(),
							rand));
		return proposalDistributions.get(rand.nextInt(proposalDistributions
				.size()));
	}
	
	private EvolutionParameterProposalDistribution nextEvolParaProposal(Random rand) {
		if (evolProposalDistributions.isEmpty())
			evolProposalDistributions.addAll(EvolutionParameterProposalDistribution.Util
					.proposalList(evolProposalOptions, rand));
		return evolProposalDistributions.get(rand.nextInt(evolProposalDistributions
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