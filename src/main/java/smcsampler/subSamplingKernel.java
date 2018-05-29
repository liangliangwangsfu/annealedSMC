package smcsampler;
import static nuts.util.CollUtils.list;

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

import cern.colt.Arrays;
import ev.ex.DataGenerator;
import ev.ex.TreeGenerators;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.mcmc.ProposalDistribution;
//import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import pty.smc.models.CTMC;
import pty.smc.models.FastDiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import fig.prob.Distrib;
import fig.prob.Gamma;
import goblin.Taxon;
import ma.MSAParser;
import ma.MSAPoset;
import ma.SequenceType;
import nuts.math.Graph;
import nuts.math.HashGraph;
import nuts.util.CollUtils;
import nuts.util.Indexer;
import pepper.Encodings;

public class subSamplingKernel implements subSMCsamplerKernel<UnrootedTreeState>{
	@Option public static double AnnealDeltaProposalRate = 10.0;    
	private int nAnnealing = 500;
	@Option public static boolean printBranchLengthMagnitudes = false;
	@Option private Dataset data = null;
	private final UnrootedTreeState initial;
	private double temperature = 0;
	private double newtemperature = 0;
	private double defaultTemperatureDifference = 0;
	private boolean initializing = true;
	private int currentIter=0;
	private int nTemperatures = 0;
	private Pair<List<Integer>, Pair<Integer, Double>> temperatureIndex = new Pair<List<Integer>, Pair<Integer, Double>>(null, null);
	//private Pair<List<Integer>, Pair<Integer, Double>> temperatureDiffIndex = new Pair<List<Integer>, Pair<Integer, Double>>(null, null);
	private Map<Taxon, LikelihoodModelCalculator> leaves4one = null;
	private Map<Taxon, LikelihoodModelCalculator> leaves4digits = null;
	private int tempdigitsIndex = 0;
	private Map<Taxon, LikelihoodModelCalculator> leaves4oneDiff = null;
	private Map<Taxon, LikelihoodModelCalculator> leaves4digitsDiff4Data1 = null;
	private Map<Taxon, LikelihoodModelCalculator> leaves4digitsDiff4Data2 = null;
	private int tempdigitsIndexDiff = 0;
	private List<Integer> temperatureDiffIndex4One = null;
	private Pair<List<Integer>, List<Double>> temperatureDiffIndex4Digits = null;

	private LinkedList<ProposalDistribution> proposalDistributions = null;
	private ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;

	public double getTemperature() {
		return temperature;
	}
	
//	public void setTemperature(double temperature) {
//		this.temperature = temperature;
//	}

	public void setTemperatureIndex(Pair<List<Integer>, Pair<Integer, Double>> temperatureIndex) {
		this.temperatureIndex = temperatureIndex;
	}
	
	public void setTemperatureDiffIndex4One(List<Integer> temperatureDiffIndex4One) {
		this.temperatureDiffIndex4One = temperatureDiffIndex4One;
	}
	
	public void setTemperatureDiffIndex4Digits(Pair<List<Integer>, List<Double>> temperatureDiffIndex4Digits) {
		this.temperatureDiffIndex4Digits = temperatureDiffIndex4Digits;
	}
	
	  public void setleaves4digitsDiff4Data1(Map<Taxon, LikelihoodModelCalculator> leaves4digitsDiff4Data1) {
		  this.leaves4digitsDiff4Data1 = leaves4digitsDiff4Data1;
	  }
	  public void setleaves4digitsDiff4Data2(Map<Taxon, LikelihoodModelCalculator> leaves4digitsDiff4Data2) {
		  this.leaves4digitsDiff4Data2 = leaves4digitsDiff4Data2;
	  }
	  public void setleaves4oneDiff(Map<Taxon, LikelihoodModelCalculator> leaves4oneDiff) {
		  this.leaves4oneDiff = leaves4oneDiff;
	  }
	  public void setleaves4digits(Map<Taxon, LikelihoodModelCalculator> leaves4digits) {
		  this.leaves4digits = leaves4digits;
	  }
	  public void setleaves4one(Map<Taxon, LikelihoodModelCalculator> leaves4one) {
		  this.leaves4one = leaves4one;
	  }
	  
	  public void setTempdigitsIndexDiff(int tempdigitsIndexDiff) {
		  this.tempdigitsIndexDiff = tempdigitsIndexDiff;
	  }
	  
	  public void setTempdigitsIndex(int tempdigitsIndex) {
		  this.tempdigitsIndex = tempdigitsIndex;
	  }
	
	public void setData(Dataset data) {
		this.data = data;
	}

	private double temperatureDifference = 0;
	

	public double getTemperatureDifference() {
		return temperatureDifference;
	}
//	public void setTemperatureDifference(double temperatureDifference) {				
//		double newTemp=temperature + temperatureDifference;
//		if(newTemp>=1.0) 
//		{
//			newtemperature=1.0;
//			this.temperatureDifference=1.0-temperature;			
//		}
//		else
//		{
//			newtemperature=newTemp;
//			this.temperatureDifference=temperatureDifference;
//		}			
//		//		LogInfo.logs("temperature: " + temperature +"  newtemperature: " + newtemperature);
//	}
	
	Map<Taxon, LikelihoodModelCalculator> Templeaves4one(Dataset data, List<Integer> temperature4OneIndex){
		int sitesLength = temperature4OneIndex.size();
		List<Taxon> leafNames = new ArrayList<Taxon>();
		Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();
		Map<Taxon, double[][]> observations = data.observations();
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(sitesLength, 2);
		for (Taxon lang : observations.keySet()) {
			double[][] observationSite = new double[sitesLength][data.nCharacter(1)];
			for(int i = 0; i < sitesLength ; i++)
				observationSite[i] = observations.get(lang)[temperature4OneIndex.get(i)];
			
			leafNames.add(lang);
			leaves.put(
			lang,
			FastDiscreteModelCalculator.observation(ctmc,
					observationSite, false));
		}
		return(leaves);
	}
	
	Map<Taxon, LikelihoodModelCalculator> Templeaves4digits(Dataset data, int temperature4DigitsIndex){
		List<Taxon> leafNames = new ArrayList<Taxon>();
		Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();
		Map<Taxon, double[][]> observations = data.observations();
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(1, 2);
		for (Taxon lang : observations.keySet()) {
			double[][] observationSite = new double[1][data.nCharacter(1)];
			observationSite[0] = observations.get(lang)[temperature4DigitsIndex];
			
			leafNames.add(lang);
			leaves.put(
			lang,
			FastDiscreteModelCalculator.observation(ctmc,
					observationSite, false));
		}
		return(leaves);
	}
	
	public subSamplingKernel(UnrootedTreeState initial, int nTemperatures, LinkedList<ProposalDistribution> proposalDistributions, ProposalDistribution.Options proposalOptions) 
	{
		this.initial = initial;		 
		this.proposalDistributions=proposalDistributions;
		this.proposalOptions=proposalOptions;
		this.nTemperatures=nTemperatures;
	}
	public UnrootedTreeState getInitial() { return initial; }

	public boolean isLastIter()
	{
		return temperatureIndex.getFirst().size() > nTemperatures;
	}

	public UnrootedTreeState sampleFromPrior(Random rand, UnrootedTreeState current)
	{ 
		Gamma exponentialPrior = Gamma.exponential(AnnealDeltaProposalRate);
		//		RootedTree proprosedRTree = TreeGenerators.sampleExpNonclock(rand, current.getUnrootedTree().leaves(), AnnealDeltaProposalRate);
		//		UnrootedTreeState proposedState = current.copyAndChange(UnrootedTree.fromRooted(proprosedRTree));
		//		RootedTree proprosedRTree = generate(rand, exponentialPrior,current.getUnrootedTree().leaves()); 		
		UnrootedTreeState proposedState = current.copyAndChange(generate(rand, exponentialPrior, current.getUnrootedTree().leaves()));				
		return proposedState;
	}

	@Override
	public Pair<UnrootedTreeState, Double> next(
			Random rand,
			UnrootedTreeState current)
	{
		if (currentIter==0) 
		{
			current = sampleFromPrior(rand, current);
			return Pair.makePair(current, 0.0); 
		}
		//tempdigitsIndex = temperatureIndex.getFirst().get(temperatureIndex.getFirst().size())+1;
		
		ProposalDistribution proposal = nextProposal(rand);         
		UnrootedTree currenturt=current.getNonClockTree();
		Pair<UnrootedTree,Double> result = proposal.propose(currenturt, rand);
		UnrootedTreeState  proposedState = null;
		//computation of log(w) is different ...//
		//double logw = temperatureDifference * current.logLikelihood();
		double logw = 0.0;
		//compute the sites for weights
		//tempdigitsIndexDiff = temperatureDiffIndex.getSecond().getFirst();
		
		if(temperatureDiffIndex4Digits.getFirst().size() > 0 && temperatureDiffIndex4One.size() == 0) {
			logw = UnrootedTreeState.computeLogLikelihood(currenturt,leaves4digitsDiff4Data1)*temperatureDiffIndex4Digits.getSecond().get(0);
			if(temperatureDiffIndex4Digits.getFirst().size() > 1 && temperatureDiffIndex4One.size() == 0) {
				logw = logw + UnrootedTreeState.computeLogLikelihood(currenturt,leaves4digitsDiff4Data2)*temperatureDiffIndex4Digits.getSecond().get(1);
			}
		}
		if(temperatureDiffIndex4Digits.getFirst().size() == 0 && temperatureDiffIndex4One.size() > 0) {
		    logw = UnrootedTreeState.computeLogLikelihood(currenturt,leaves4oneDiff);
		}
		if(temperatureDiffIndex4Digits.getFirst().size() > 0 && temperatureDiffIndex4One.size() > 0) {
		    logw = UnrootedTreeState.computeLogLikelihood(currenturt,leaves4oneDiff) + UnrootedTreeState.computeLogLikelihood(currenturt,leaves4digitsDiff4Data1)*temperatureDiffIndex4Digits.getSecond().get(0);
			if(temperatureDiffIndex4Digits.getFirst().size() > 1 && temperatureDiffIndex4One.size() > 0) {
				logw = logw + UnrootedTreeState.computeLogLikelihood(currenturt,leaves4digitsDiff4Data2)*temperatureDiffIndex4Digits.getSecond().get(1);
			}
		}
		
		
		if (result != null) // might happen e.g. when trying to do nni with 3 leaves
		{
			double logTargetDenCurrent= 0.0;
			double logLikRatio = 0.0;
			proposedState = current.copyAndChange(result.getFirst());
			final double logProposalRatio = result.getSecond();
			tempdigitsIndex = temperatureIndex.getSecond().getFirst();
			
			//three cases: all one leaves, just digits, one leaves plus digits...
			if(tempdigitsIndex == -1) {
				//leaves4one = Templeaves4one(data, temperatureIndex.getFirst());
				logTargetDenCurrent = logTargetDensityOne(leaves4one, current);				
				logLikRatio = logTargetDensityOne(leaves4one, proposedState) -logTargetDenCurrent;  	
			}
			if(temperatureIndex.getFirst().size() == 0) {
				//leaves4digits = Templeaves4digits(data, tempdigitsIndex);
				logTargetDenCurrent = logTargetDensitydigits(leaves4digits, temperatureIndex.getSecond().getSecond(), current);				
				logLikRatio = logTargetDensitydigits(leaves4digits, temperatureIndex.getSecond().getSecond(), proposedState) -logTargetDenCurrent;  	
			}
			if((temperatureIndex.getFirst().size() > 0)&&(tempdigitsIndex >= 0)) {
				//leaves4one = Templeaves4one(data, temperatureIndex.getFirst());
				//leaves4digits = Templeaves4digits(data, tempdigitsIndex);
				logTargetDenCurrent = logTargetDensityOneDigits(leaves4one, leaves4digits, temperatureIndex.getSecond().getSecond(), current);				
				logLikRatio = logTargetDensityOneDigits(leaves4one, leaves4digits, temperatureIndex.getSecond().getSecond(), proposedState) -logTargetDenCurrent;  				
			}
						
			final double ratio = Math.min(1,
					Math.exp(logProposalRatio + logLikRatio));
			if (Double.isNaN(ratio))
				throw new RuntimeException();
			if (rand.nextDouble() >= ratio) {
				proposedState = current;
			}
		}
		temperature = newtemperature;
		return Pair.makePair(proposedState, logw);    
	}

	public double logTargetDensityOne(Map<Taxon, LikelihoodModelCalculator> leaves4one, UnrootedTreeState uts)
	{
		return uts.getLogPrior()+ UnrootedTreeState.computeLogLikelihood(uts.getNonClockTree(),leaves4one);	
	}
	
	public double logTargetDensityOneDigits(Map<Taxon, LikelihoodModelCalculator> leaves4one, Map<Taxon, LikelihoodModelCalculator> leaves4digits, double temp4digits, UnrootedTreeState uts)
	{
		return uts.getLogPrior()+ UnrootedTreeState.computeLogLikelihood(uts.getNonClockTree(),leaves4one) + temp4digits*UnrootedTreeState.computeLogLikelihood(uts.getNonClockTree(),leaves4digits);	
	}
	
	public double logTargetDensitydigits(Map<Taxon, LikelihoodModelCalculator> leaves4digits, double temp4digits, UnrootedTreeState uts)
	{
		return uts.getLogPrior() + temp4digits*UnrootedTreeState.computeLogLikelihood(uts.getNonClockTree(),leaves4digits);		
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
		
		
		//List<Taxon> leaves = MSAParser.parseMSA(alignmentInputFile).taxa();
		
		//List<Taxon> leaves = msa.taxa();
		//Taxon taxa = leaves.get(8);
		System.out.println(data.nSites());
		System.out.println(msa);
		
		String msaString = msa.toString();
		
		int Ntemp = 500;
		double[] originalTemp = new double[Ntemp];
		for(int i = 0; i < Ntemp; i++)
			originalTemp[i] = Math.pow((i+1)/Ntemp, 3);
		//Map<Taxon, double[][]> observations = data.observations();
		//char[] column1 = msaString.charAt(1,1);
		//double[][] dataMatrix = observations.get();
//		double[][] dataMat = observations.get(taxa);
//		List<Taxon> list = new ArrayList<Taxon>(observations.keySet());
//		double [][] cObs = observations.get(taxa);
		int[] sites = new int[len-1];
		for(int i = 0; i < len-1; i++)
			sites[i] = i;
		double[] temp = {1, 1, 0.1};
		List<Taxon> leafNames = new ArrayList<Taxon>();
		Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();
		Map<Taxon, double[][]> observations = data.observations();
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(sites.length, 2);
		for (Taxon lang : observations.keySet()) {
			double[][] observationSite = new double[sites.length][data.nCharacter(1)];
			for(int i = 0; i < sites.length ; i++)
				observationSite[i] = observations.get(lang)[sites[i]];
			leafNames.add(lang);
			//System.out.println(Arrays.toString(observationSite));
//			leaves.put(
//					lang,
//					FastDiscreteModelCalculator.observation(ctmc,
//							observations.get(lang), false));
			leaves.put(
			lang,
			FastDiscreteModelCalculator.observation(ctmc,
					observationSite, false));
			//System.out.println(lang);
			//System.out.println(Arrays.toString(observations.get(lang)[0]));
			
		}
		
		double loglikelihood = UnrootedTreeState.computeLogLikelihood(urt,leaves);
		
		int l1 = 1001;
		int l2 = 10;
		int nIndex = 0;
		if(l1 == (l1/l2)*l2) {
			nIndex = l1/l2;
		}else {
			nIndex = l1/l2 + 1;
		}
		System.out.println(nIndex);

		



		
	}



}
