package evolmodel;

import java.io.File;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Random;
import conifer.trees.StandardNonClockPriorDensity;
import ev.ex.DataGenerator;
import fig.basic.Option;
import fig.basic.Pair;
import fig.prob.Gamma;
import goblin.Taxon;
import ma.MSAParser;
import ma.MSAPoset;
import ma.SequenceType;
import nuts.maxent.SloppyMath;
import nuts.util.Indexer;
import pepper.Encodings;
import pty.RootedTree;
import pty.io.Dataset;
import ev.ex.TreeGenerators;
import pty.UnrootedTree;
import pty.mcmc.*;
import pty.smc.models.CTMC;



public class SteppingStoneEvolPara implements Runnable{
	@Option public File alignmentInputFile = null;
	@Option public SequenceType st = SequenceType.DNA;
	@Option public int nSamplesEachChain = 4000;
	@Option public int nChains = 13;
	@Option public double alpha = 0.3;
	@Option public static double GammapriorRatio = 10.0;
	@Option public double csmc_trans2tranv=2.0;
	@Option public static EvolutionModel evolModel = EvolutionModel.K2P;
	@Option public static  LinkedList<EvolutionParameterProposalDistribution> evolProposalDistributions = null;
	@Option public static  EvolutionParameterProposalDistribution.Options evolProposalOptions = EvolutionParameterProposalDistribution.Util._defaultProposalDistributionOptions;
	@Option public static  Dataset dataset;
	@Option public Random rand = new Random(1);

//	public static class Options{
//		@Option public Prior prior = Prior.EXP;
//		@Option public Random rand = new Random(1);
//	}

	private int nSamples = (int)(nSamplesEachChain*0.8);
	private int nburn = nSamplesEachChain - nSamples;
	private double logZ = 0.0;



	public void setnSamplesEachChain(int nSamplesEachChain) {
		this.nSamplesEachChain = nSamplesEachChain;
		this.nSamples = (int)(nSamplesEachChain*0.8);
		this.nburn = nSamplesEachChain - nSamples;
	}


	public static Gamma exponentialPrior = Gamma.exponential(GammapriorRatio);
	public static StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
			exponentialPrior);
	@Override
	public void run()
	{
		SteppingStoneEvolPara newrun = new SteppingStoneEvolPara();
		int nSamples=(int) (nSamplesEachChain*0.8);
		int nburn=nSamplesEachChain-nSamples;    
		System.out.println(nChains);
		System.out.println(nSamples);
		System.out.println(nSamplesEachChain);
		System.out.println(csmc_trans2tranv);
		MSAPoset align = MSAParser.parseMSA(alignmentInputFile);
		Dataset data = Dataset.DatasetUtils.fromAlignment(align, st);
		EvolutionParameters evolPara = null;
		if(evolModel==EvolutionModel.K2P)
			evolPara = new EvolutionParameters.K2P(csmc_trans2tranv);
		else
			evolPara = new EvolutionParameters.GTR(new double[]{0.26, 0.18, 0.17, 0.15, 0.11, 0.13, 0.25, 0.25, 0.25, 0.25});  // GTR model				
		CTMC ctmc = evolModel.instantiateCTMC(evolPara, data.nSites());	
		UnrootedTree initTree = initTree(new Random(3), align.taxa());
		UnrootedTreeState initUrootedTreeState = UnrootedTreeState.initFastState(initTree, dataset, ctmc, priorDensity);
		newrun.SteppingStone(rand, new UnrootedTreeEvolParameterState(initUrootedTreeState,evolPara), data, ctmc, nChains, nSamplesEachChain, alpha);
		double logZ = newrun.getNormalizer();
		estimateNormalizer(logZ);
		//newrun.estimateNormalizer(logZ);		 
	}


	public static PhyloSampler.Options _defaultPhyloSamplerOptions = new PhyloSampler.Options();	

	private static double logNormalizer(List<Double> Loglikelihood, double temperature1, double temperature2, int nSamples) {
		double out = 0.0;
		double[] logterm = new double[nSamples];
		for(int i = 0; i < nSamples; i++) {
			logterm[i] = (temperature2 - temperature1)*Loglikelihood.get(i);
			//System.out.println("temperature: "+ temperature1 + " Loglikelihood: "+ Loglikelihood.get(i));
		}
		out = SloppyMath.logAdd(logterm) - Math.log(1.0*nSamples);
		return(out);
	}


	private static Pair<UnrootedTree, List<Double>>propagation(Random r, double temperature, int nburn, int nSamples, CTMC ctmc, StandardNonClockPriorDensity priorDensity, UnrootedTreeEvolParameterState initState) {
		List<UnrootedTree> proposedSample = new ArrayList();
		List<Double> proposedLoglikelihood = new ArrayList();
	//	Random r = new Random();
		//we use the tree from previous chain //
		/*		UnrootedTree sample = initTree(r, msa.taxa());
		UnrootedTreeState temp = UnrootedTreeState.initFastState(sample, data, ctmc, priorDensity);*/
		UnrootedTreeEvolParameterState temp = initState;
		UnrootedTreeEvolParameterState proposedState = null;

		ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
		List<ProposalDistribution> proposalDistributions = ProposalDistribution.Util.proposalList(proposalOptions, initState.getUnrootedTreeState().getNonClockTree(), r);

		EvolutionParameterProposalDistribution.Options evolProposalOptions = EvolutionParameterProposalDistribution.Util._defaultProposalDistributionOptions;
		if(evolModel==EvolutionModel.GTR) {
			evolProposalOptions.useGTRProposal = true;
			evolProposalOptions.useK2PProposal = false;
		}		
		List<EvolutionParameterProposalDistribution> evolProposalDistributions  = EvolutionParameterProposalDistribution.Util.proposalList(evolProposalOptions, r);

		
		for(int i = 0; i < (nburn + nSamples); i++) {
			ProposalDistribution nextProposal = proposalDistributions.get(r.nextInt(proposalDistributions.size()));
			EvolutionParameterProposalDistribution nextEvolParaProposal = evolProposalDistributions.get(r.nextInt(evolProposalDistributions.size()));  	
			proposedState = proposal(r, temp, nextProposal, nextEvolParaProposal, temperature);
			temp = proposedState;
			if(i >= nburn){
				proposedLoglikelihood.add(proposedState.getUnrootedTreeState().getLogLikelihood());		
			}
		}
		return Pair.makePair(proposedState.getUnrootedTreeState().getNonClockTree(), proposedLoglikelihood);
	}

	public static  UnrootedTree initTree(Random rand, List<Taxon> leaves){	
		return UnrootedTree.fromRooted(TreeGenerators.sampleExpNonclock(rand,leaves, 10.0));				
	}


	private static EvolutionParameterProposalDistribution nextEvolParaProposal(Random rand) {
		if (evolProposalDistributions.isEmpty())
			evolProposalDistributions.addAll(EvolutionParameterProposalDistribution.Util
					.proposalList(evolProposalOptions, rand));
		return evolProposalDistributions.get(rand.nextInt(evolProposalDistributions
				.size()));
	}


	private static double logTargetDensity(double temperature, UnrootedTreeState uts)
	{
		return uts.getLogPrior()+ temperature*uts.getLogLikelihood();	
	}

	public static UnrootedTreeEvolParameterState proposal(Random rand, UnrootedTreeEvolParameterState current, ProposalDistribution proposal, EvolutionParameterProposalDistribution evolParaProposal, double temperature){
		UnrootedTreeEvolParameterState proposedState = null;
		UnrootedTreeState proposedUrootedTreeState = null;    
		UnrootedTree currenturt=current.getUnrootedTreeState().getNonClockTree();
		Pair<UnrootedTree,Double> result = proposal.propose(currenturt, rand);
		if (result != null) // might happen e.g. when trying to do nni with 3 leaves
		{
			double logTargetDenCurrent = logTargetDensity(temperature, current.getUnrootedTreeState());
//			EvolutionParameterProposalDistribution evolParaProposal = nextEvolParaProposal(rand);  
			Pair<EvolutionParameters, Double> evolProposalRe = evolParaProposal.propose(current.getEvolParameter(), rand);
			EvolutionParameters proposedEvolPar = evolProposalRe.getFirst();    // proposed evolutionary parameters
			CTMC ctmc = evolModel.instantiateCTMC(proposedEvolPar, dataset.nSites());			
			proposedUrootedTreeState = UnrootedTreeState.initFastState(result.getFirst(), dataset, ctmc, priorDensity);  // result.getFirst() gives the proposed tree and ctmc is updated with the newly proposed evolutionary parameters			
			final double logProposalRatio = result.getSecond()+evolProposalRe.getSecond();		
			double logLikRatio = logTargetDensity(temperature, proposedUrootedTreeState) - logTargetDenCurrent;  
			final double ratio = Math.min(1,
					Math.exp(logProposalRatio + logLikRatio));
			if (Double.isNaN(ratio))
				throw new RuntimeException();
			if (rand.nextDouble() >= ratio) {
				proposedState = current;			
			}else
			{
				proposedState = new UnrootedTreeEvolParameterState(proposedUrootedTreeState, proposedEvolPar);
			}
		}
		return proposedState;
	}

	/*	public static enum SStempScheme{
		Beta{
			public double[] generateTemp(int N, double alpha) {
				double[] Temp = new double[N];
				for(int i = 0; i < N; i++) {
					Temp[i] = Math.pow(i*1.0/((N-1)*1.0), 1.0/alpha);
				}
				return(Temp);				
			}
		},
		Evenly{
			public double[] generateTemp(int N,  double alpha){
				double[] Temp = new double[N];
				for(int i = 0; i < N; i++) {
					Temp[i] = i*1.0/((N-1)*1.0);
				}
				return(Temp);		
			}
		};

		public abstract double[] generateTemp(int N, double alpha); 

		public double[] SSTemp(int N, double alpha){

			return generateTemp(N, alpha);
		}

	}	*/

	public static enum SStempScheme{
		Beta{
			public double[] generateTemp(int N, double alpha) {
				double[] Temp = new double[N];
				for(int i = 0; i < N; i++) {
					Temp[i] = Math.pow(i*1.0/((N-1)*1.0), 1.0/alpha);
				}
				return(Temp);				
			}
		},
		Evenly{
			public double[] generateTemp(int N,  double alpha){
				double[] Temp = new double[N];
				for(int i = 0; i < N; i++) {
					Temp[i] = i*1.0/((N-1)*1.0);
				}
				return(Temp);		
			}
		};

		public abstract double[] generateTemp(int N, double alpha); 

		public double[] SSTemp(int N, double alpha){

			return generateTemp(N, alpha);
		}

	}	

	public static double max(double[] vec){
		int N = vec.length;
		double max = -999999;
		for(int i = 0; i < N; i++) {
			if(vec[i] > max) {
				max = vec[i];
			}
		}
		return max;		
	}


	public UnrootedTreeEvolParameterState sampleFromPrior(Random rand, UnrootedTreeEvolParameterState current)
	{ 
		Gamma exponentialPrior = Gamma.exponential(GammapriorRatio);
		UnrootedTreeState proposedState = current.getUnrootedTreeState().copyAndChange(AnnealingKernelTreeEvolPara.generate(rand, exponentialPrior, current.getUnrootedTreeState().getUnrootedTree().leaves()));
		current.getEvolParameter().sampleFromPrior(rand);		
		return new UnrootedTreeEvolParameterState(proposedState, current.getEvolParameter());
	}


	public void SteppingStone(Random r, UnrootedTreeEvolParameterState current,  Dataset data, CTMC ctmc, int nChains, int nSamplesEachChain, double alpha) {
		setnSamplesEachChain(nSamplesEachChain);
		Pair<UnrootedTree, List<Double>> proposedState = null;
		List<Double> proposedLoglikelihood = new ArrayList<Double>();
		final double[] temperatureSchedule = SStempScheme.Beta.generateTemp(nChains, alpha);
		//		UnrootedTree initTree = null;
		//		UnrootedTreeState initTreeState = null;
		UnrootedTreeEvolParameterState initState = null;  
		for(int i = 0; i < nSamples; i++) {
		//	Random r =  new Random();
			initState = sampleFromPrior(r, current);
			proposedLoglikelihood.add(initState.getUnrootedTreeState().getLogLikelihood());
		}

		logZ = logZ + logNormalizer(proposedLoglikelihood, temperatureSchedule[0], temperatureSchedule[1], nSamples);

		System.out.println(logZ);

		for(int t = 2; t < nChains; t++) {
			proposedState = propagation(r, temperatureSchedule[t-1], nburn, nSamples, ctmc, priorDensity, initState);
			logZ = logZ + logNormalizer(proposedState.getSecond(), temperatureSchedule[t-1], temperatureSchedule[t], nSamples);
			//			initTree = proposedState.getFirst();
			//			initTreeState =  UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);
			System.out.println(logZ);
		}
	}

	public void setnChains(int nChains) {
		this.nChains = nChains;
	}

	public double getNormalizer() {
		return logZ;
	}

	public void estimateNormalizer(double logZ) {
		this.logZ = logZ;
	}


	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int nTaxa = 10;
		double treeRate = 10;
		int len = 100;
		boolean generateDNAdata = true;
		Random randTree =  new Random(1);
		Random randData =  new Random(2);

		SequenceType sequenceType = SequenceType.DNA;
		RootedTree rt = TreeGenerators.sampleExpNonclock(randTree, nTaxa, treeRate);
		UnrootedTree urt = UnrootedTree.fromRooted(rt);
		System.out.println(urt);

		Indexer<Character> indexer; 
		indexer = Encodings.dnaEncodings().nonGapCharactersIndexer();	
		DataGenerator generator= new DataGenerator(generateDNAdata, rt);  

		MSAPoset  msa = generator.generate(randData, len);
		Dataset data = Dataset.DatasetUtils.fromAlignment(msa, sequenceType);
		//	CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(data.nSites(), 2);

		SteppingStoneEvolPara newrun = new SteppingStoneEvolPara();
		dataset = data; 
	    evolModel = EvolutionModel.K2P;
		int nChains = 50;
		int nSamplesEachChain  = 1000;
		double alpha = 1.0/3.0;


		EvolutionParameters evolPara = null;
		if(evolModel==EvolutionModel.K2P)
			evolPara = new EvolutionParameters.K2P(2.0);
		else
			evolPara = new EvolutionParameters.GTR(new double[]{0.26, 0.18, 0.17, 0.15, 0.11, 0.13, 0.25, 0.25, 0.25, 0.25});  // GTR model				
		CTMC ctmc = evolModel.instantiateCTMC(evolPara, data.nSites());	
		UnrootedTree initTree = initTree(new Random(3), msa.taxa());
		UnrootedTreeState initUrootedTreeState = UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);		
		newrun.SteppingStone(new Random(), new UnrootedTreeEvolParameterState(initUrootedTreeState,evolPara), data, ctmc, nChains, nSamplesEachChain, alpha);
		double logZ = newrun.getNormalizer();
		System.out.println(logZ);



	}



}
