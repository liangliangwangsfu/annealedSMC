package smcsampler;

import java.io.File;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Random;

import bayonet.distributions.Multinomial;
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
import pty.RandomRootedTrees;
import pty.RootedTree;
import pty.RootedTree.RootingInfo;
import pty.io.Dataset;
import pty.io.Dataset.DatasetUtils;
import ev.ex.TreeGenerators;
import ev.ex.PhyloSamplerMain.SStempScheme;
import ev.poi.processors.TreeDistancesProcessor;
import pty.UnrootedTree;
import pty.mcmc.*;
import pty.mcmc.PhyloSampler.Prior;
import pty.smc.models.CTMC;
import pty.smc.models.CTMC.GTRIGammaCTMC;



public class LinkedImportanceSampling implements Runnable{
	@Option public File alignmentInputFile = null;
	@Option public SequenceType st = SequenceType.DNA;
	@Option public int nChains = 13;
	@Option public double alpha = 0.3;
	@Option public double GammapriorRatio = 10.0;
	@Option public int nSamplesEachChain = 4000;
	@Option public double csmc_trans2tranv=2.0;
	private int nSamples = (int)(nSamplesEachChain*0.8);
	private int nburn = nSamplesEachChain - nSamples;
	
	public void setnSamplesEachChain(int nSamplesEachChain) {
		this.nSamplesEachChain = nSamplesEachChain;
		this.nSamples = (int)(nSamplesEachChain*0.8);
		this.nburn = nSamplesEachChain - nSamples;
	}
	
	public static class Options{
	    @Option public Prior prior = Prior.EXP;
	    @Option public Random rand = new Random(1);
	}
		
	
	private double logZ = 0.0;

	private Gamma exponentialPrior = Gamma.exponential(GammapriorRatio);
	private StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
			exponentialPrior);
	 @Override
	  public void run()
	  {
		 LinkedImportanceSampling newrun = new LinkedImportanceSampling();
		 int nSamples=(int) (nSamplesEachChain*0.8);
		 int nburn=nSamplesEachChain-nSamples;    
		 System.out.println(nChains);
		 System.out.println(nSamples);
		 System.out.println(csmc_trans2tranv);
		 MSAPoset align = MSAParser.parseMSA(alignmentInputFile);
		 Dataset data = Dataset.DatasetUtils.fromAlignment(align, st);
		 CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(data.nSites(), csmc_trans2tranv);	
		 newrun.LinkedIS(align, data, ctmc, nChains, nSamplesEachChain, alpha);
		 double logZ = newrun.getNormalizer();
		 estimateNormalizer(logZ);
		 //newrun.estimateNormalizer(logZ);		 
	  }

	
	public static PhyloSampler.Options _defaultPhyloSamplerOptions = new PhyloSampler.Options();
	
	
	private static int LinkedIndex(List<Double> sampleLoglikelihood, double temperature1, double temperature2, int nSamples) {
		int index = -1;
		double[] probability = new double[nSamples];
		double[] logtemp = new double[nSamples];
		Random r = new Random();
		for(int i = 0; i < nSamples; i++) {
			logtemp[i] = sampleLoglikelihood.get(i)*(temperature2 - temperature1)/2;
		}
		probability = normalizeWeights(logtemp);
		index = Multinomial.sampleMultinomial(r, probability); 
		return(index);
	}
	
	
	private static double logNormalizer(List<Double> Loglikelihood1, List<Double> Loglikelihood2, double temperature1, double temperature2, int nSamples) {
		double out = 0.0;
		double[] logterm1 = new double[nSamples];
		double[] logterm2 = new double[nSamples];
		for(int i = 0; i < nSamples; i++) {
			logterm1[i] = (temperature2 - temperature1)/2*Loglikelihood1.get(i);
			logterm2[i] = (temperature1 - temperature2)/2*Loglikelihood2.get(i);		
		}	
		out = SloppyMath.logAdd(logterm1) - SloppyMath.logAdd(logterm2);
		return(out);
	}
	
	
	private static Pair<List<UnrootedTree>, List<Double>>propagation(UnrootedTree sample, double sampleLoglikelihood, double temperature, int nburn, int nSamples, CTMC ctmc, Dataset data, StandardNonClockPriorDensity priorDensity) {
		List<UnrootedTree> proposedSample = new ArrayList();
		List<Double> proposedLoglikelihood = new ArrayList();
		//proposedLoglikelihood.add(sampleLoglikelihood);
		Random r = new Random();
		UnrootedTreeState temp = UnrootedTreeState.initFastState(sample, data, ctmc, priorDensity);

		UnrootedTreeState proposedState = null;
		
		ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
		//List<ProposalDistribution> proposalDistributions = ProposalDistribution.Util.proposalList(proposalOptions, r, temperature);
		List<ProposalDistribution> proposalDistributions = ProposalDistribution.Util.proposalList(proposalOptions, temp.getNonClockTree(), r);
		
		for(int i = 0; i < (nburn + nSamples); i++) {
			ProposalDistribution nextProposal = proposalDistributions.get(r.nextInt(proposalDistributions
					.size()));
			//System.out.println(nextProposal);
			proposedState = proposal(r, temp, nextProposal, temperature);
			temp = proposedState;
			if(i >= nburn - 1){
				proposedSample.add(proposedState.getNonClockTree());
				proposedLoglikelihood.add(proposedState.getLogLikelihood());		
			}
				
		}
		return Pair.makePair(proposedSample, proposedLoglikelihood);
		
	}
	
	public static  UnrootedTree initTree(Random rand, List<Taxon> leaves){	
		return UnrootedTree.fromRooted(TreeGenerators.sampleExpNonclock(rand,leaves, 10.0));				
	}
	

	
	private static double logTargetDensity(double temperature, UnrootedTreeState uts)
	{
		return uts.getLogPrior()+ temperature*uts.getLogLikelihood();	
	}
	
	public static UnrootedTreeState proposal(Random rand, UnrootedTreeState current, ProposalDistribution proposal, double temperature){
		UnrootedTreeState proposedState = null;
		double logPosterior = 0.0;   
		UnrootedTree currenturt=current.getNonClockTree();
		Pair<UnrootedTree,Double> result = proposal.propose(currenturt, rand);
		if (result != null) // might happen e.g. when trying to do nni with 3 leaves
		{
			//double logTargetDenCurrent=temperature*current.logLikelihood();
			double logTargetDenCurrent= logTargetDensity(temperature, current);
			proposedState = current.copyAndChange(result.getFirst());
			final double logProposalRatio = result.getSecond();
			//double logLikRatio = temperature*proposedState.getLogLikelihood() - logTargetDenCurrent;  
			double logLikRatio = logTargetDensity(temperature, proposedState) - logTargetDenCurrent;  
			//System.out.println("temperature = "+ temperature);
			//System.out.println("LogLikRatio = "+ logLikRatio + "; logProposalRatio = "+ logProposalRatio);
			final double ratio = Math.min(1,
					Math.exp(logProposalRatio + logLikRatio));
			if (Double.isNaN(ratio))
				throw new RuntimeException();
			if (rand.nextDouble() >= ratio) {
				proposedState = current;
			}
		}
		return proposedState;
	}
	
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
	
	
	public static double[] normalizeWeights(double[] logweights) {
		int N = logweights.length;
		double[] normalizedweights = new double[N];
		double maxLogweights  = max(logweights);
		double[] centralizeWeights = new double[N];
		double csumweights = 0;
		for(int i = 0; i < N; i++) {
			centralizeWeights[i] = Math.exp(logweights[i] - maxLogweights);
			csumweights = csumweights + centralizeWeights[i];
		}
		
		for(int i = 0; i < N; i++) {
			normalizedweights[i] = centralizeWeights[i]/csumweights;
		}
		
		return normalizedweights;
	}
	

	public void LinkedIS(MSAPoset  msa, Dataset data, CTMC ctmc, int nChains, int nSamplesEachChain, double alpha) {
		setnSamplesEachChain(nSamplesEachChain);
		Random r =  new Random();
		Pair<List<UnrootedTree>, List<Double>> proposedState = null;
		Pair<List<UnrootedTree>, List<Double>> proposedState1 = null;
		List<UnrootedTree> proposedSample = new ArrayList<UnrootedTree>();
		List<UnrootedTree> proposedSample1 = new ArrayList<UnrootedTree>();
		List<Double> proposedLoglikelihood = new ArrayList<Double>();
		List<Double> proposedLoglikelihood1 = new ArrayList<Double>();

		final double[] temperatureSchedule = SStempScheme.Beta.generateTemp(nChains, alpha);
		UnrootedTree LinkedTree = initTree(r, msa.taxa());
		UnrootedTree initTree = null;
		UnrootedTreeState initTreeState = null;
		double LinkedTreeLoglikelihood = 0.0;
		int I = (int)(r.nextDouble()*nSamples);

		System.out.println(nSamplesEachChain);
		System.out.println(nburn);
		System.out.println(nSamples);
		for(int i = 0; i < nSamples; i++) {
			initTree = initTree(r, msa.taxa());
			proposedSample.add(LinkedTree);	
		    initTreeState = UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);
			proposedLoglikelihood.add(initTreeState.getLogLikelihood());
		}
		LinkedTree = proposedSample.get(I);
		UnrootedTreeState LinkedState = UnrootedTreeState.initFastState(LinkedTree, data, ctmc, priorDensity);
		
		double t1 = temperatureSchedule[1];
		proposedState1 = propagation(LinkedTree, LinkedState.getLogLikelihood(), t1, nburn, nSamples, ctmc, data, priorDensity);
		
		logZ = logZ + logNormalizer(proposedLoglikelihood, proposedState1.getSecond(), temperatureSchedule[0], temperatureSchedule[1], nSamples);
		
		System.out.println(logZ);

		for(int t = 2; t < nChains; t++) {
			proposedState = proposedState1;
			t1 = temperatureSchedule[t];
			I = LinkedIndex(proposedState.getSecond(), temperatureSchedule[t-1], t1, nSamples);
			LinkedTree = proposedState.getFirst().get(I);
			LinkedTreeLoglikelihood = proposedState.getSecond().get(I);
			proposedState1 = propagation(LinkedTree, LinkedTreeLoglikelihood, t1, nburn, nSamples, ctmc, data, priorDensity);
			logZ = logZ + logNormalizer(proposedState.getSecond(), proposedState1.getSecond(), temperatureSchedule[t-1], t1, nSamples);
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
		int nTaxa = 30;
		double treeRate = 10;
		int len = 200;
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
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(data.nSites(), 2);
		
/*		for(int i = 0; i < T; i++ ) {
			System.out.println(temperatureSchedule[i]);
		}*/
		LinkedImportanceSampling newrun = new LinkedImportanceSampling();
		int nChains = 50;
		int nSamplesEachChain  = 1000;
		double alpha = 1.0/3.0;
		newrun.LinkedIS(msa, data, ctmc, nChains, nSamplesEachChain, alpha);
		double logZ = newrun.getNormalizer();
		System.out.println(logZ);


		
	}



}
