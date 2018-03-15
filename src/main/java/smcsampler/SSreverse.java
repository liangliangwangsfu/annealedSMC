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


public class SSreverse implements Runnable{
	@Option public File alignmentInputFile = null;
	@Option public SequenceType st = SequenceType.DNA;
	@Option public int nSamplesEachChain = 4000;
	@Option public int nChains = 13;
	@Option public double alpha = 0.3;
	@Option public double GammapriorRatio = 10.0;
	
	
	public static class Options{
	    @Option public Prior prior = Prior.EXP;
	    @Option public Random rand = new Random(1);
	}
		
	private int nSamples = (int)(nSamplesEachChain*0.8);
	private int nburn = nSamplesEachChain - nSamples;
	private double logZ = 0.0;
	
	public void setnSamplesEachChain(int nSamplesEachChain) {
		this.nSamplesEachChain = nSamplesEachChain;
		this.nSamples = (int)(nSamplesEachChain*0.8);
		this.nburn = nSamplesEachChain - nSamples;
	}
	

	private Gamma exponentialPrior = Gamma.exponential(GammapriorRatio);
	private StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
			exponentialPrior);
	 @Override
	  public void run()
	  {
		 SteppingStone newrun = new SteppingStone();
		 int nSamples=(int) (nSamplesEachChain*0.8);
		 int nburn=nSamplesEachChain-nSamples;    
		 System.out.println(nChains);
		 System.out.println(nSamples);
		 System.out.println(nSamplesEachChain);
		 MSAPoset align = MSAParser.parseMSA(alignmentInputFile);
		 Dataset data = Dataset.DatasetUtils.fromAlignment(align, st);
		 CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(data.nSites(), 1);	
		 newrun.SteppingStone(align, data, ctmc, nChains, nSamplesEachChain, alpha);
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
	
	//private static Pair<List<UnrootedTree>, List<Double>>propagation(MSAPoset  msa, double temperature, int nburn, int nSamples, CTMC ctmc, Dataset data, StandardNonClockPriorDensity priorDensity, UnrootedTreeState initTree) {
	private static Pair<UnrootedTree, List<Double>>propagation(MSAPoset  msa, double temperature, int nburn, int nSamples, CTMC ctmc, Dataset data, StandardNonClockPriorDensity priorDensity, UnrootedTreeState initTree) {
		List<UnrootedTree> proposedSample = new ArrayList();
		List<Double> proposedLoglikelihood = new ArrayList();
		Random r = new Random();
		//we use the tree from previous chain //
/*		UnrootedTree sample = initTree(r, msa.taxa());
		UnrootedTreeState temp = UnrootedTreeState.initFastState(sample, data, ctmc, priorDensity);*/
		UnrootedTreeState temp = initTree;
		UnrootedTreeState proposedState = null;
		
		ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
		List<ProposalDistribution> proposalDistributions = ProposalDistribution.Util.proposalList(proposalOptions, initTree.getNonClockTree(), r);
		
		for(int i = 0; i < (nburn + nSamples); i++) {
			ProposalDistribution nextProposal = proposalDistributions.get(r.nextInt(proposalDistributions
					.size()));
			proposedState = proposal(r, temp, nextProposal, temperature);
			temp = proposedState;
			if(i >= nburn){
				proposedLoglikelihood.add(proposedState.getLogLikelihood());		
			}
				
		}
		return Pair.makePair(proposedState.getNonClockTree(), proposedLoglikelihood);
		
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
			double logTargetDenCurrent = logTargetDensity(temperature, current);
			proposedState = current.copyAndChange(result.getFirst());
			final double logProposalRatio = result.getSecond();
			//double logLikRatio = temperature*proposedState.getLogLikelihood() - logTargetDenCurrent;  
			double logLikRatio = logTargetDensity(temperature, proposedState) - logTargetDenCurrent;  
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
					Temp[i] = Math.pow(i*1.0/(N*1.0), 1.0/alpha);
				}
				return(Temp);				
			}
		},
		Evenly{
			public double[] generateTemp(int N,  double alpha){
				double[] Temp = new double[N];
				for(int i = 0; i < N; i++) {
					Temp[i] = i*1.0/(N*1.0);
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
	

	public void SteppingStone(MSAPoset  msa, Dataset data, CTMC ctmc, int nChains, int nSamplesEachChain, double alpha) {
		setnSamplesEachChain(nSamplesEachChain);
		//Random r =  new Random();
		//Pair<List<UnrootedTree>, List<Double>> proposedState = null;
		Pair<UnrootedTree, List<Double>> proposedState = null;
		//List<UnrootedTree> proposedSample = new ArrayList<UnrootedTree>();
		List<Double> proposedLoglikelihood = new ArrayList<Double>();

		final double[] temperatureSchedule = SStempScheme.Beta.generateTemp(nChains, alpha);
		UnrootedTree initTree = null;
		UnrootedTreeState initTreeState = null;
		
		//First sample from the posterior chain...
		Random r0 =  new Random();
		initTree = initTree(r0, msa.taxa());
		initTreeState = UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);
		proposedState = null;
		proposedState = propagation(msa, 1.0, nburn, nSamples, ctmc, data, priorDensity, initTreeState);
		initTree = proposedState.getFirst();
		initTreeState = UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);
		
		for(int t = 0; t < nChains-1; t++) {
			proposedState = null;
			//t1 = temperatureSchedule[t-1];
			proposedState = propagation(msa, temperatureSchedule[nChains-1-t], nburn, nSamples, ctmc, data, priorDensity, initTreeState);
			if(t == 0) {
				logZ = logZ + logNormalizer(proposedState.getSecond(), temperatureSchedule[nChains-1-t], 1, nSamples);
			}
			else {
				logZ = logZ + logNormalizer(proposedState.getSecond(), temperatureSchedule[nChains-1-t], temperatureSchedule[nChains-t], nSamples);
			}
			
			//initTree = proposedState.getFirst().get(nSamples - 1);
			initTree = proposedState.getFirst();
			initTreeState =  UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);
			System.out.println(logZ);
		}


		for(int i = 0; i < nSamples; i++) {
			Random r =  new Random();
			initTree = initTree(r, msa.taxa());
		    initTreeState = UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);
			proposedLoglikelihood.add(initTreeState.getLogLikelihood());
		}

		logZ = logZ + logNormalizer(proposedLoglikelihood, temperatureSchedule[0], temperatureSchedule[1], nSamples);
		
		System.out.println(logZ);
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
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(data.nSites(), 2);
		
		SSreverse newrun = new SSreverse();
		int nChains = 50;
		int nSamplesEachChain  = 2000;
		double alpha = 1.0/3.0;
		newrun.SteppingStone(msa, data, ctmc, nChains, nSamplesEachChain, alpha);
		double logZ = newrun.getNormalizer();
		System.out.println(logZ);


		
	}



}

