package pty.smc;
import java.util.LinkedList;
import java.util.Random;

import nuts.math.Sampling;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.mcmc.NonClockTreeState;
import pty.mcmc.PhyloSampler;
import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import ev.ex.TreeGenerators;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import goblin.Taxon;


public class AnnealingKernel implements ParticleKernel<UnrootedTreeState>
{
	@Option public static double AnnealDeltaProposalRate = 10.0;   //  
	// @Option public static boolean useOptimal = true;
	// @Option public static boolean useLazy = false;

	@Option
	public static int nAnnealing = 500;

	@Option
	public static boolean useSteppingStone = true;
	
	@Option public static boolean printBranchLengthMagnitudes = false;
	private final UnrootedTreeState initial;
	// private int nIter=1000;
	private boolean lastIter = false;

	private double temperature = 0;
	private double newtemperature = 0;

	public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	private double temperatureDifference = 1.0 / nAnnealing;

	public double getTemperatureDifference() {
		return temperatureDifference;
	}
	public void setTemperatureDifference(double temperatureDifference) {
		this.temperatureDifference = temperatureDifference;
		newtemperature = Math.min(temperature + temperatureDifference,
				1.0);
		// LogInfo.logsForce("temperature " + temperature
		// + " temperatureDifference " + temperatureDifference + " "
		// + "newtemperature: " + newtemperature);
		LogInfo.logs("newtemperature: " + newtemperature);

		if (newtemperature == 1.0)
			lastIter = true;

	}
	private int currentIter=0;
	private LinkedList<ProposalDistribution> proposalDistributions = null;
	private ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
	public AnnealingKernel(UnrootedTreeState initial, LinkedList<ProposalDistribution> proposalDistributions, ProposalDistribution.Options proposalOptions) 
	{
		this.initial = initial;		 
		this.proposalDistributions=proposalDistributions;
		this.proposalOptions=proposalOptions; 
	}
	public UnrootedTreeState getInitial() { return initial; }


	
	public boolean isLastIter()
	{
		return lastIter;
	}
	
	public Object _next(
			Random rand,
			UnrootedTreeState current, boolean isPeek)
 {
		if (temperature == 0) {
			RootedTree proprosedRTree = TreeGenerators.sampleCoalescent(rand,
					current.getUnrootedTree().leaves(), false);
			UnrootedTreeState proposedState = current
					.copyAndChange(UnrootedTree.fromRooted(proprosedRTree));
			temperature = newtemperature;
			Double logw = temperatureDifference
					* (proposedState.logLikelihood() + prior(proposedState
							.getNonClockTree()));
			return Pair.makePair(proposedState, logw);
		}
		//
		PhyloSampler  sampler = new PhyloSampler();
		// sampler.setOutputText(true);
		sampler.init(current);
		// sampler.sample();
		// double oldtemperature = (1.0*(currentIter-1))/nAnnealing, temperature
		// = (1.0*currentIter)/nAnnealing;
		// System.out.println("temperature: "+temperature+"; currentIter"+temperature*nAnnealing);
		// temperatureDifference = temperature - oldtemperature;
		// System.out.println("newtemperature " + temperatureDifference);
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
				// final double energyLogRatio = sampler.energy(proposedState) -
				// sampler.currentEnergy();
				// final double ratio = Math.min(1,Math.exp(logProposalRatio -
				// energyLogRatio/temperature));

				double currentPrior = current.getLogPrior(), proposePrior = proposedState
						.getLogPrior();
						
						//prior(current.getUnrootedTree()), proposePrior = prior(proposedState
						//.getUnrootedTree()); // log

				// System.out.println("currentPrior " + currentPrior + " "
				// + current.getLogPrior());
				//
				// System.out.println("proposePrior " + proposePrior + " "
				// + proposedState.getLogPrior());

				double logLikRatio = newtemperature
						* (proposePrior + proposedState.getLogLikelihood())
						// - temperature
						- newtemperature
						* (currentPrior + current.getLogLikelihood());
				final double ratio = Math.min(1,
						Math.exp(logProposalRatio + logLikRatio));

				if (Double.isNaN(ratio))
					throw new RuntimeException();
				if (rand.nextDouble() >= ratio) {
					proposedState = current;
					// System.out.println("ratio: "+ratio );
				}
			}
		}
		temperature = newtemperature;
		// Double logw =
		// (1.0/(nIter-currentIter+1)-1.0/(nIter-currentIter+2))*(current.logLikelihood());
		// //+prior(current.getNonClockTree());
	 double logw =0; 
	 if(useSteppingStone) logw = temperatureDifference * current.logLikelihood();
	 else	 
			logw = temperatureDifference
					* (current.logLikelihood() + prior(current
							.getNonClockTree()));

		//System.out.println("nIter: "+nIter+" currentIter:" + currentIter + (nIter-currentIter)+":  logw:"+(1.0/(nIter-currentIter+1)-1.0/(nIter-currentIter+2)));

		// Math.log(0.5*(lambda*Math.exp(-lambda*delta)+0.1))
		return Pair.makePair(proposedState, logw);    
	}



	@Override
	public Pair<UnrootedTreeState, Double> next(Random rand,
			UnrootedTreeState current)
			{
		return (Pair) _next(rand, current, false);
			}

	@Override
	public int nIterationsLeft(UnrootedTreeState partialState) {	
		// System.out.println("nIterationsLeft: "+ (nIter-currentIter));
		// return (nIter-currentIter);
		return (nAnnealing);
	}

	public void setCurrentIter(int  currentIter) {	
		this.currentIter=currentIter; 
	}

	public int getCurrentIter() {
		return this.currentIter;
	}

	private void sample(Random rand, NonClockTreeState currentState) {
	}

	private ProposalDistribution nextProposal(Random rand) {
		if (proposalDistributions.isEmpty())
			proposalDistributions.addAll(ProposalDistribution.Util
					.proposalList(proposalOptions, initial.getNonClockTree(),
							rand));
		return proposalDistributions.get(rand.nextInt(proposalDistributions
				.size()));
	}

	private double prior(UnrootedTree urt) {
		double result = 0.0;
		for (UnorderedPair<Taxon, Taxon> edge : urt.edges()) {
			// System.out.println(edge
			// + " "
			// + urt.branchLength(edge)
			// + "---"
			// + Sampling.exponentialLogDensity(AnnealDeltaProposalRate,
			// urt.branchLength(edge)));

			// result += Sampling.exponentialLogDensity(AnnealDeltaProposalRate,
			// urt.branchLength(edge)); // TODO: check the other places where I
			// used this function. I should use 1.0 / AnnealDeltaProposalRate
			// rather than AnnealDeltaProposalRate.
			result += Sampling.exponentialLogDensity(
					1.0 / AnnealDeltaProposalRate, urt.branchLength(edge));
		}
		return result;
	}

}