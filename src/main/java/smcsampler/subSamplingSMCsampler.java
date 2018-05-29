package smcsampler;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;

import conifer.trees.StandardNonClockPriorDensity;
import ev.ex.DataGenerator;
import ev.ex.TreeGenerators;
import ev.poi.processors.TreeDistancesProcessor;
import monaco.process.ProcessSchedule;
import nuts.io.CSV;
import nuts.lang.ArrayUtils;
import nuts.math.Fct;
import nuts.math.Id;
import nuts.math.MeasureZeroException;
import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.Hasher;
import nuts.util.Indexer;
import pepper.Encodings;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import pty.smc.PartialCoalescentState;
import pty.smc.models.CTMC;
import pty.smc.models.FastDiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import smcsampler.SMCSampler.AdaptiveScheme;
import smcsampler.SMCSampler.ParticleMapperProcessor;
import smcsampler.SMCSampler.ParticleProcessor;
import smcsampler.SMCSampler.ResamplingStrategy;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.Parallelizer;
import fig.prob.Gamma;
import fig.prob.Multinomial;
import goblin.BayesRiskMinimizer;
import goblin.Taxon;
import goblin.BayesRiskMinimizer.LossFct;
import ma.MSAParser;
import ma.MSAPoset;
import ma.SequenceType;

public final class subSamplingSMCsampler<S> {
	@Option
	public boolean verbose = false;
	public int N = 100;
	@Option
	public Random rand = new Random(1);
	@Option
	public boolean resampleLastRound = false;
	@Option
	public int nThreads = 1;
	@Option
	public ResamplingStrategy resamplingStrategy = ResamplingStrategy.ESS;
	@Option
	public int nSites = 100;
	@Option
	public int nFakeSites = 100;
	@Option
	public int nSitesPerIndex = 10;
	@Option
	public Dataset data = null;

	@Option
	public  AdaptiveScheme adaptiveScheme=AdaptiveScheme.CONSTANT; 
	public static double essRatioThreshold = 0.5;

	private ProcessSchedule schedule = null;
	private double ess = N;
	
	private List<Integer> tempdigitsIndexDiff;
	
	private int tempdigitsIndex;
	//private List<Double> tempDiffList = new ArrayList<Double>();
	//	private double cess=1.0;
	private double tempDiff=0;
	
	private List<Double> DeterministictemperatureDifference = null;
	
	private List<Double> Deterministictemperature = null;
	
	//This is where we revised
	//private List<double[]> newTemperature = null;
	
	//private List<double[]> newTemperatureDiff = null;
	
	private double[][] newTemperature;
	
	private double[][] newTemperatureDiff;
	
	
//	public List<Double> getDeterministicTemperature(){
//		return tempDiffList;
//	}
	
	private Pair<List<Integer>, Pair<Integer, Double>> temperatureIndex = null;
	
	private List<Integer> temperatureDiffIndex4One = null;
	
	private Pair<List<Integer>, List<Double>> temperatureDiffIndex4Digits = null;
	
	private Map<Taxon, LikelihoodModelCalculator> leaves4one = null;
	private Map<Taxon, LikelihoodModelCalculator> leaves4digits = null;
	private Map<Taxon, LikelihoodModelCalculator> leaves4oneDiff = null;
	private Map<Taxon, LikelihoodModelCalculator> leaves4digitsDiff4Data1 = null;
	private Map<Taxon, LikelihoodModelCalculator> leaves4digitsDiff4Data2 = null;

	
	public void setData(Dataset data) {
		this.data = data;
	}
	
	public void setnSites(int nSites) {
		this.nSites = nSites;
	}
	
	public void setnSitesPerIndex(int nSitesPerIndex){
		this.nSitesPerIndex = nSitesPerIndex;
	}
	
	public void setnFakeSites(int nFakeSites) {
		this.nFakeSites = nFakeSites;
	}
	
	public void setDeterministicTemperatureDifference(List<Double> temperatureDifference) {
		this.DeterministictemperatureDifference = temperatureDifference;
	}
	
	public void setDeterministicTemperature(List<Double> temperature) {
		this.Deterministictemperature = temperature;
	}
	
	
//	public void setNewTemperature(List<Double> temperature, int nSites){
//		int ntemperatures = temperature.size();
//		List<double[]> Temperatures = new ArrayList<double[]>();
//		for(int i = 0; i < ntemperatures; i++) {
//			double[] TempVec = new double[nSites];
//			for(int j = 0; j < nSites; j++) {
//				TempVec[j] = psiFunction(j, nSites, temperature.get(i));
//			}
//			Temperatures.add(TempVec);
//		}
//		
//		List<double[]> TemperatureDiff = new ArrayList<double[]>();
//		TemperatureDiff.add(Temperatures.get(0));
//		for(int i = 1; i < Temperatures.size(); i++) {
//			System.out.println(i);
//			double[] TempVecDiff = new double[nSites];
//			for(int j = 0; j < nSites; j++) {
//				TempVecDiff[j] = Temperatures.get(i)[j] - Temperatures.get(i-1)[j];
//			}
//			TemperatureDiff.add(TempVecDiff);
//		}
//		
//		this.newTemperature = Temperatures;
//		this.newTemperatureDiff = TemperatureDiff;
//	}
	
	public void setNewTemperature2(List<Double> temperature, int nSites){
		int ntemperatures = temperature.size();
		double[][] Temperatures = new double[ntemperatures][nSites];
		for(int i = 0; i < ntemperatures; i++) {
			for(int j = 0; j < nSites; j++) {
				Temperatures[i][j] = psiFunction(j, nSites, temperature.get(i));
			}
		}

		this.newTemperature = Temperatures;
	}
	
	private double[] NewTemperatureDiff2(int i, int nSites){
		double[] TemperatureDiff = new double[nSites];
		if(i == 0) {
			TemperatureDiff = newTemperature[0];
		}else {
			for(int j = 0; j < nSites; j++) {
				//TempVecDiff[j] = Temperatures.get(i)[j] - Temperatures.get(i-1)[j];
				TemperatureDiff[j] = newTemperature[i][j] - newTemperature[i-1][j];
			}
		}
		//this.newTemperature = Temperatures;
		return(TemperatureDiff);
	}
	
	public void setNewTemperatureDiff(int nSites){
		int ntemperatures = newTemperature.length;

		double[][] TemperatureDiff = new double[ntemperatures][nSites];
		TemperatureDiff[0] = newTemperature[0];
		System.out.println(ntemperatures);
		System.out.println(nSites);

		for(int i = 1; i < ntemperatures; i++) {
			System.out.println(i);
			//double[] TempVecDiff = new double[nSites];
			for(int j = 0; j < nSites; j++) {
				//TempVecDiff[j] = Temperatures.get(i)[j] - Temperatures.get(i-1)[j];
				TemperatureDiff[i][j] = newTemperature[i][j] - newTemperature[i-1][j];
			}
			//TemperatureDiff.add(TempVecDiff);
		}
		
		//this.newTemperature = Temperatures;
		this.newTemperatureDiff = TemperatureDiff;
	}
	
//	private void setNewTemperatureDiff(List<double[]> Temperatures, int nSites){
//		List<double[]> TemperatureDiff = new ArrayList<double[]>();
//		TemperatureDiff.add(Temperatures.get(0));
//		for(int i = 1; i < Temperatures.size(); i++) {
//			double[] TempVecDiff = new double[nSites];
//			for(int j = 0; j < nSites; j++) {
//				TempVecDiff[j] = Temperatures.get(i)[j] - Temperatures.get(i-1)[j];
//			}
//			TemperatureDiff.add(TempVecDiff);
//		}
//		this.newTemperatureDiff = TemperatureDiff;
//	}
	
	public double[] getOneTemperatureSeq(int i){
		//return newTemperature.get(i);
		return newTemperature[i];
	}
	
	public double[] getOneTemperatureDiffSeq(int i){
		//return newTemperatureDiff.get(i);
		return newTemperatureDiff[i];
	}
	
	private double psiFunction(int ss, int nsites, double phi){
		double rt = 0.0;
		int s = ss + 1;
		if(phi >= 1.0*s/(1.0*nsites)){
			rt = 1;
		}
		if(phi <= 1.0*(s-1)/(1.0*nsites)){
			rt = 0;
		}
		if((phi < 1.0*s/(1.0*nsites))&&(phi > 1.0*(s-1)/(1.0*nsites))){
			rt = phi*nsites - (s-1);
		}
		return(rt);
	}
/*	public double getDeterministicTemperatureDifference(int t) {
		return(DeterministictemperatureDifference.get(t));
	}*/
	
	private boolean useCESS=true;

	public boolean isUseCESS() {
		return useCESS;
	}

	public void setUseCESS(boolean useCESS) {
		this.useCESS = useCESS;
	}

	@Option
	public boolean adaptiveTempDiff = false;
	public static double alpha = 0.90;

	public int adaptiveType = 0;
	public PrintWriter smcSamplerOut = null;
	public PrintWriter smcSamplerOut2 = null;

	public double getEss() {
		return ess;
	}

	public void setEss(double ess) {
		this.ess = ess;
	}

	public void setProcessSchedule(ProcessSchedule schedule) {
		this.schedule = schedule;
	}
	
	List<Integer> temperatureIndexReal(List<Integer> temperature4OneIndex, int sitesLength, int nSitesPerIndex){
		List<Integer> temperatureIndexReal = new ArrayList<Integer>();
		temperatureIndexReal.add(temperature4OneIndex.get(0)*nSitesPerIndex);
		for(int i = 1; i < sitesLength; i++) {
			temperatureIndexReal.add(temperatureIndexReal.get(i-1)+1);
		}
		return(temperatureIndexReal);
	}
	
	List<Integer> temperatureIndexDigitsReal(int temperature4DigitsIndex, int sitesLength, int nSitesPerIndex){
		List<Integer> temperatureIndexReal = new ArrayList<Integer>();
		temperatureIndexReal.add(temperature4DigitsIndex*nSitesPerIndex);
		for(int i = 1; i < sitesLength; i++) {
			temperatureIndexReal.add(temperatureIndexReal.get(i-1)+1);
		}
		return(temperatureIndexReal);
	}
	
	Map<Taxon, LikelihoodModelCalculator> Templeaves4one(Dataset data, List<Integer> temperature4OneIndex){
		int sitesLength = 1;
		int fakelength = temperature4OneIndex.size();
		List<Integer> temperature4OneIndexReal = null;
		if(temperature4OneIndex.get(fakelength-1) == nFakeSites) {
			sitesLength = (fakelength-1)*nSitesPerIndex + nSites - (nFakeSites-1)*nSitesPerIndex;
		}else {
			sitesLength = fakelength*nSitesPerIndex;			
		}
		temperature4OneIndexReal = temperatureIndexReal(temperature4OneIndex, sitesLength, nSitesPerIndex);
		//int sitesLength = temperature4OneIndex.size();
		List<Taxon> leafNames = new ArrayList<Taxon>();
		Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();
		Map<Taxon, double[][]> observations = data.observations();
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(sitesLength, 2);
		for (Taxon lang : observations.keySet()) {
			double[][] observationSite = new double[sitesLength][data.nCharacter(1)];
			for(int i = 0; i < sitesLength ; i++)
				observationSite[i] = observations.get(lang)[temperature4OneIndexReal.get(i)];
				//observationSite[i] = observations.get(lang)[temperature4OneIndex.get(i)];
			
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
		int sitesLength = 1;
		if(temperature4DigitsIndex == nFakeSites) {
			sitesLength = nSites - (nFakeSites-1)*nSitesPerIndex;
		}else {
			sitesLength = nSitesPerIndex;			
		}
		List<Integer> temperature4DigitsIndexReal = null;
		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(sitesLength, 2);
		temperature4DigitsIndexReal = temperatureIndexDigitsReal(temperature4DigitsIndex, sitesLength, nSitesPerIndex);
		for (Taxon lang : observations.keySet()) {
			double[][] observationSite = new double[sitesLength][data.nCharacter(1)];
			for(int i = 0; i < sitesLength ; i++)
				observationSite[i] = observations.get(lang)[temperature4DigitsIndexReal.get(i)];
			
			//double[][] observationSite = new double[1][data.nCharacter(1)];
			//observationSite[0] = observations.get(lang)[temperature4DigitsIndex];
			
			leafNames.add(lang);
			leaves.put(
			lang,
			FastDiscreteModelCalculator.observation(ctmc,
					observationSite, false));
		}
		return(leaves);
	}
	
	Pair<List<Integer>, Pair<Integer, Double>> tempSeqProcess(double[] temperature){
		List<Integer> Left = new ArrayList<Integer>();
		Integer RightInt = -1;
		Double RightDigit = 0.0;
		//Pair<Integer, Double> Right = null;
//		if(temperature[0] < 1) {
//			Right = temperature[0];
//		}else {
			for(int i = 0; i < temperature.length; i++) {
				if(temperature[i] == 1.0) {Left.add(i);}
				if(temperature[i] < 1.0 && temperature[i] >0.0) {RightInt = i; RightDigit = temperature[i];}
			}
		//}	
		Pair<Integer, Double> Right = new Pair <Integer, Double>(RightInt,RightDigit);		
		return(new Pair <List<Integer>, Pair<Integer, Double>>(Left,Right) );
	}
	
	List<Integer> tempSeqProcessDiff4One(double[] temperatureDiff){
		List<Integer> Left = new ArrayList<Integer>();
			for(int i = 0; i < temperatureDiff.length; i++) {
				if(temperatureDiff[i] == 1.0) {Left.add(i);}
			}

		return(Left);
	}
	
	Pair<List<Integer>, List<Double>> tempSeqProcessDiff4digits(double[] temperatureDiff){
		List<Integer> Left = new ArrayList<Integer>();
		List<Double> Right = new ArrayList<Double>();
			for(int i = 0; i < temperatureDiff.length; i++) {
				if(temperatureDiff[i] < 1.0 && temperatureDiff[i] >0.0) {Left.add(i); Right.add(temperatureDiff[i]);}
			}

		return(new Pair<>(Left, Right));
	}





	public static enum ResamplingStrategy {
		ALWAYS {
			@Override
			public boolean needResample(double[] w) {
				return true;
			}
		},
		NEVER {
			@Override
			public boolean needResample(double[] w) {
				return false;
			}
		},
		ESS {
			@Override
			public boolean needResample(final double[] weights) {
				final double ess = ess(weights);
				final double threshold = essRatioThreshold * weights.length;
				return ess < threshold;
			}
		};
		abstract boolean needResample(double[] weights);
	}

	public static enum AdaptiveScheme {
		CONSTANT {
			@Override
			double alpha(double currentAlpha, int t){
				return alpha;
			}
		},
		Adaptive1 {
			@Override
			double alpha(double currentAlpha,int t){
				return currentAlpha + 0.5 * Math.pow(1 - alpha, 2.0)
				* Math.pow(alpha, t);
			}
		},
		Adaptive2 {
			@Override
			double alpha(double currentAlpha,int t){
				return currentAlpha + 0.5 * Math.pow(1 - alpha, 3.0)
				* (t + 1) * Math.pow(alpha, t);
			}
		};				
		abstract double alpha(double currentAlpha, int t);
	}

	public static double ess(double[] ws) {
		double sumOfSqr = 0.0;
		for (double w : ws)sumOfSqr += w * w;
		return 1.0/sumOfSqr;
	}

	public static double cess(double[] normalizedPreviousWs, double[] logws) {
		if(normalizedPreviousWs.length!=logws.length)  
			throw new RuntimeException("normalizedPreviousWs and ws should have the same length!");
		double [] numerator=new double[logws.length];
		double [] denominator =new double[logws.length];
		for (int i=0;i<logws.length;i++){
			denominator[i] = Math.log(normalizedPreviousWs[i])+2*logws[i];
			numerator[i] = Math.log(normalizedPreviousWs[i])+logws[i];
		}		
		double tmp=SloppyMath.logAdd(numerator);
		return Math.exp(2*tmp-SloppyMath.logAdd(denominator));
	}

	public static <S> void bootstrapFilter(final subSMCsamplerKernel<S> kernel, 
			final TreeDistancesProcessor processor, final int N,
			final Random rand) {
		subSamplingSMCsampler<S> PF = new subSamplingSMCsampler<S>();
		PF.N = N;
		PF.rand = rand;
		PF.resampleLastRound = false;
		PF.resamplingStrategy = ResamplingStrategy.ESS;
		PF.sample(kernel, processor);
	}

	private long[] seeds;
	private List<S> samples;
	private double[] logWeights;
	private double[] normalizedWeights;
	private double[] logIncrementalWeights;
	private double varLogZ = 0;
	public List<S> getSamples() {
		return samples;
	}
	public double[] getLogWeights() {
		return logWeights;
	}
	
	public double[] getLogIncrementalWeights() {
		return logIncrementalWeights;
	}
	
	private void propagateAndComputeWeights(final subSMCsamplerKernel<S> kernel,
			final int t)
	{		
		if (verbose)
			LogInfo.track("Processing...", false);
		seeds = Sampling.createSeeds(N, rand); // so that result for a given
		// randomization is
		// #-of-threads-invariant
		Parallelizer<Integer> parallelizer = new Parallelizer<Integer>(nThreads);
		parallelizer.setPrimaryThread();
		parallelizer.process(CollUtils.ints(N),
				new Parallelizer.Processor<Integer>() {
			public void process(Integer x, int _i, int _n, boolean log) {
				if (log && ((x + 1) % 5 == 0))
					if (verbose)
						LogInfo.logs("Particle " + (x + 1) + "/" + N);
				Random rand = new Random(seeds[x]);
				final Pair<S, Double> current = kernel.next(rand,
						samples.get(x));
				if (current == null) {
					samples.set(x, null);
					logWeights[x] = Double.NEGATIVE_INFINITY;
				} else {
					samples.set(x, current.getFirst());
					//logIncrementalWeights[x] = current.getSecond();
					logWeights[x] = Math.log(normalizedWeights[x])+ current.getSecond();				
				}
			}
		});		
		if(verbose) LogInfo.end_track();		
	}
	private double lognorm = 0.0;
	public double estimateNormalizer() {
		return lognorm;
	}

	public double estimateNormalizerVariance() {
		return varLogZ;
	}

	private void init(final subSMCsamplerKernel<S> kernel)
	{
		lognorm = 0.0;
		kernel.setData(data);
		samples = new ArrayList<S>(N);
		S initial = kernel.getInitial();
		for (int n = 0; n < N; n++) samples.add(initial);
		logWeights = new double[N]; // init to zero	
		logIncrementalWeights = new double[N]; // init to zero	
		normalizedWeights = new double[N];
		for(int i=0;i<normalizedWeights.length;i++) normalizedWeights[i]=1.0/N;
	}

//	public double temperatureDifference(final double alpha,double absoluteAccuracy,double min,double max) {
//		final double[] logLike = new double[samples.size()];
//		for (int n = 0; n < samples.size(); n++) {
//			UnrootedTreeState urt = ((UnrootedTreeState) samples.get(n));
//			logLike[n] = urt.getLogLikelihood(); 
//		}
//		int maxEval = 100;
//		UnivariateFunction f = new UnivariateFunction() {
//			public double value(double x) {
//				double[] logWeightLikePriorVec = new double[samples.size()];
//				if(useCESS)
//				{   
//					for (int n = 0; n < samples.size(); n++) 	logWeightLikePriorVec[n] = logLike[n]* x; //incremental weights										
//					return (cess(normalizedWeights,logWeightLikePriorVec) - alpha);
//				}
//				else
//				{					
//					for (int n = 0; n < samples.size(); n++) 	logWeightLikePriorVec[n] = Math.log(normalizedWeights[n])+logLike[n]* x; 	
//					NumUtils.expNormalize(logWeightLikePriorVec); // no need to normalize 
//					return (ess(logWeightLikePriorVec) - alpha*ess);
//				}
//
//			}
//		};
//		double result = 0;
//		try {
//			final double relativeAccuracy = absoluteAccuracy * 0.0001;
//			// final double absoluteAccuracy0 = 1.0e-6;
//			PegasusSolver solver = new PegasusSolver(relativeAccuracy,
//					absoluteAccuracy);
//			result = solver.solve(maxEval, (UnivariateFunction) f, min, max);
//			//			LogInfo.logsForce("Solver successful!");
//		} catch (RuntimeException e) {
//			LogInfo.logsForce("Solver Fail!");
//			result = -1;
//		}			
//		return result;
//	}


	/**
	 * @param <S>
	 * @param kernel
	 *            model
	 * @param tdp
	 *            what to do with the produced sample
	 */
	public void sample(final subSMCsamplerKernel<S> kernel,
			final TreeDistancesProcessor tdp) {
		init(kernel);
		if(smcSamplerOut!=null)smcSamplerOut.println(CSV.header("t", "ESS", "tempDiff", "temp"));
		//if(smcSamplerOut!=null)smcSamplerOut.println(CSV.header("tempDiff"));
		double alpha0 = alpha;
		//double[] newTemperatureDiff4T = new double[nSites];
		double[] newTemperatureDiff4T = new double[nFakeSites];
		//	for (int t = 0; t < T && !isLastIter; t++) {
		int t=0;
//		temperatureIndex = tempSeqProcess(newTemperature.get(t));
//		temperatureDiffIndex = tempSeqProcess(newTemperatureDiff.get(t));
//		kernel.setTemperatureIndex(temperatureIndex);
//		kernel.setTemperatureDiffIndex(temperatureDiffIndex);
		//while(!kernel.isLastIter()){	
		//while(t < newTemperature.size()){	
		while(t < newTemperature.length){
			kernel.setCurrentIter(t);
			//newTemperatureDiff4T = NewTemperatureDiff2(t, nSites);
			newTemperatureDiff4T = NewTemperatureDiff2(t, nFakeSites);
//			if (adaptiveTempDiff) {
//				alpha0 =adaptiveScheme.alpha(alpha0, t); 
//				if (t > 0) {						
//					tempDiff = temperatureDifference(alpha, 1.0e-10, 0, 0.2);
//					if(tempDiff == -1) tempDiff = kernel.getDefaultTemperatureDifference();
//					tempDiffList.add(tempDiff);
//				}
//			} else {
//				//tempDiff = kernel.getDefaultTemperatureDifference();  
//				tempDiff = DeterministictemperatureDifference.get(t); 
//				//System.out.println(tempDiff);
//			}
			//temperatureIndex = tempSeqProcess(newTemperature.get(t));
			//temperatureDiffIndex4One = tempSeqProcessDiff4One(newTemperatureDiff.get(t));
			//temperatureDiffIndex4Digits = tempSeqProcessDiff4digits(newTemperatureDiff.get(t));
			temperatureIndex = tempSeqProcess(newTemperature[t]);
			temperatureDiffIndex4One = tempSeqProcessDiff4One(newTemperatureDiff4T);
			temperatureDiffIndex4Digits = tempSeqProcessDiff4digits(newTemperatureDiff4T);
			//temperatureDiffIndex = tempSeqProcessDiff(newTemperatureDiff.get(t));
//			System.out.println("iteration: "+ t);
//			System.out.println("Temperature index integer: " + temperatureIndex.getFirst());
//			System.out.println("Temperature index digits: " + temperatureIndex.getSecond());
//			System.out.println("TemperatureDiff index integer: " + temperatureDiffIndex4One);
//			System.out.println("TemperatureDiff index digits: " + temperatureDiffIndex4Digits);
			kernel.setTemperatureIndex(temperatureIndex);
			kernel.setTemperatureDiffIndex4One(temperatureDiffIndex4One);
			kernel.setTemperatureDiffIndex4Digits(temperatureDiffIndex4Digits);
			tempdigitsIndexDiff = temperatureDiffIndex4Digits.getFirst();
			//kernel.setTempdigitsIndexDiff(tempdigitsIndexDiff);
			

			if(tempdigitsIndexDiff.size() == 0) {
				leaves4oneDiff = Templeaves4one(data, temperatureDiffIndex4One);
				kernel.setleaves4oneDiff(leaves4oneDiff);
			}
			if(tempdigitsIndexDiff.size() > 0 && temperatureDiffIndex4One.size() == 0) {
				leaves4digitsDiff4Data1 = Templeaves4digits(data, tempdigitsIndexDiff.get(0));
				kernel.setleaves4digitsDiff4Data1(leaves4digitsDiff4Data1);
				if(tempdigitsIndexDiff.size() > 1 && temperatureDiffIndex4One.size() == 0) {
					leaves4digitsDiff4Data2 = Templeaves4digits(data, tempdigitsIndexDiff.get(1));
					kernel.setleaves4digitsDiff4Data2(leaves4digitsDiff4Data2);
				}
			}
			if(tempdigitsIndexDiff.size() > 0 && temperatureDiffIndex4One.size() > 0) {
				leaves4oneDiff = Templeaves4one(data, temperatureDiffIndex4One);
				kernel.setleaves4oneDiff(leaves4oneDiff);
				leaves4digitsDiff4Data1 = Templeaves4digits(data, tempdigitsIndexDiff.get(0));
				kernel.setleaves4digitsDiff4Data1(leaves4digitsDiff4Data1);
				if(tempdigitsIndexDiff.size() > 1 && temperatureDiffIndex4One.size() > 0) {
					leaves4digitsDiff4Data2 = Templeaves4digits(data, tempdigitsIndexDiff.get(1));
					kernel.setleaves4digitsDiff4Data2(leaves4digitsDiff4Data2);
				}
			}

			
			tempdigitsIndex = temperatureIndex.getSecond().getFirst();
			
			//three cases: all one leaves, just digits, one leaves plus digits...
			if(tempdigitsIndex == -1) {
				leaves4one = Templeaves4one(data, temperatureIndex.getFirst());
				kernel.setleaves4one(leaves4one);	 	
			}
			if(temperatureIndex.getFirst().size() == 0) {
				leaves4digits = Templeaves4digits(data, tempdigitsIndex);
				kernel.setleaves4digits(leaves4digits);	
			}
			if((temperatureIndex.getFirst().size() > 0)&&(tempdigitsIndex >= 0)) {
				leaves4one = Templeaves4one(data, temperatureIndex.getFirst());
				leaves4digits = Templeaves4digits(data, tempdigitsIndex);
				kernel.setleaves4one(leaves4one);
				kernel.setleaves4digits(leaves4digits);			
			}
			//kernel.setTemperatureIndex(tempSeqProcess(newTemperature.get(t)));
			//kernel.setTemperatureDiffIndex(tempSeqProcess(newTemperatureDiff.get(t)));
			//kernel.setTemperatureDifference(tempDiff);
			if (verbose)
				LogInfo.track("Particle generation " + (t + 1) , true); 
			propagateAndComputeWeights(kernel, t);		
			if(t>0) lognorm += SloppyMath.logAdd(logWeights); 					
							//lognorm += SloppyMath.logAdd(logIncrementalWeights)-Math.log(1.0*logIncrementalWeights.length);			
			normalizedWeights = logWeights.clone();
			NumUtils.expNormalize(normalizedWeights);
			if (verbose)
				LogInfo.logs("LargestNormalizedWeights="
						+ ArrayUtils.max(normalizedWeights));
			ess = ess(normalizedWeights);	
			//System.out.println(ess);
			//System.out.println(lognorm);
			//			cess = cess(normalizedWeights, incrementalLogWeights);
			//			System.out.println(cess);
			if(smcSamplerOut!=null)
			{
				smcSamplerOut.println(CSV.body(t, ess,kernel.getTemperatureDifference(), kernel.getTemperature()));
				//smcSamplerOut.println(CSV.body(kernel.getTemperatureDifference()));
				smcSamplerOut.flush();
			}
			if(smcSamplerOut2!=null)
			{
				//smcSamplerOut.println(CSV.body(t, ess,kernel.getTemperatureDifference(), kernel.getTemperature()));
				smcSamplerOut2.println(CSV.body(kernel.getTemperatureDifference()));
				smcSamplerOut2.flush();
			}

			if (verbose)
				LogInfo.logs("RelativeESS=" + ess
						/ normalizedWeights.length);

			if (!kernel.isLastIter()
					&& (hasNulls(samples) || resamplingStrategy
							.needResample(normalizedWeights))) {
				samples = resample(samples, normalizedWeights, rand);
				ess = N;				
				for(int i=0;i<normalizedWeights.length;i++) normalizedWeights[i]=1.0/N;
			}
			if (verbose)
				LogInfo.end_track();
			//if ((kernel.isLastIter() ) && schedule != null) {
			//System.out.println(t);
			//System.out.println(newTemperature.size());
			//if ((t == newTemperature.size()-1) && schedule != null) {
			if ((t == newTemperature.length-1) && schedule != null) {
				if (resampleLastRound) {
					// this might be useful when processing a lot of particles
					// is expensive
					Pair<List<S>, double[]> resampled = resampleAndPack(
							samples, normalizedWeights, rand);
					samples = resampled.getFirst();
					normalizedWeights = resampled.getSecond();					
				}
				if (verbose)
					LogInfo.track("Processing particles");
				for (int n = 0; n < normalizedWeights.length; n++) {
					if (verbose)
						LogInfo.logs("Particle " + (n + 1) + "/"
								+ normalizedWeights.length);
					if (samples.get(n) != null) {
						tdp.process(samples.get(n), normalizedWeights[n]);
						samples.set(n, null);
					}
				}
				if (verbose)
					LogInfo.end_track();
			}
			t++;
			//System.out.println(t);
			//System.out.println(newTemperature.size());
		}
		if(smcSamplerOut!=null)smcSamplerOut.close();
		if(smcSamplerOut2!=null)smcSamplerOut2.close();
	}



	public static double logAdd(double[] logV) {
		double max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < logV.length; i++) max=Math.max(logV[i], max); 
		if (max == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;	    
		double finalresult=0;
		for(int i=0;i<logV.length;i++)finalresult+=Math.exp(logV[i]-max);
		return  Math.log(finalresult)+max;
	}


	private boolean hasNulls(List<S> samples)
	{
		for (S item : samples)
			if (item == null)
				return true;
		return false;
	}

	private <S> Pair<List<S>, double[]> resampleAndPack(List<S> list,
			double[] w, Random rand)
	{
		if (!NumUtils.normalize(w))
			throw new MeasureZeroException();

		if (list.size() != w.length)
			throw new RuntimeException();

		Counter<Integer> packed = Sampling.efficientMultinomialSampling(rand,
				w, w.length);

		double[] resultWeight = new double[packed.size()];
		List<S> resultList = new ArrayList<S>(packed.size());

		int i = 0;
		for (int itemIdx : packed.keySet()) {
			S item = list.get(itemIdx);
			resultList.add(item);
			resultWeight[i++] = packed.getCount(itemIdx);
		}

		return Pair.makePair(resultList, resultWeight);
	}

	// TODO: can be more efficient, nlog(n) instead of n^2
	private List<S> resample(final List<S> list, final double[] w, Random rand) {
		if (!NumUtils.normalize(w))
			throw new MeasureZeroException();

		if (list.size() != w.length)
			throw new RuntimeException();

		Counter<Integer> packed = Sampling.efficientMultinomialSampling(rand,
				w, w.length);
		final List<S> result = new ArrayList<S>(list.size());
		for (int itemIdx : packed.keySet()) {
			S item = list.get(itemIdx);
			for (int cur = 0; cur < packed.getCount(itemIdx); cur++)
				result.add(item);
		}
		return result;
	}

	/**
	 * Weight is NOT in log scale, in contrast to SMCSamplerKernel.next()
	 * 
	 * @author bouchard
	 * 
	 * @param <S>
	 */
	public static interface ParticleProcessor<S> {
		public void process(S state, double weight);
	}

	public static class StoreProcessor<S> implements ParticleProcessor<S>
	{
		public List<S> particles = CollUtils.list();
		public List<Double> ws = CollUtils.list();

		@Override
		public void process(S state, double weight) {
			particles.add(state);
			ws.add(weight);
		}

		public S sample(Random rand) {
			final int idx = Sampling.sample(rand, ws);
			return particles.get(idx);
		}
	}

	public static class DoNothingProcessor<S> implements ParticleProcessor<S> {
		public void process(S state, double w) {
		}
	}




	public static class ForkedProcessor<S> implements ParticleProcessor<S>
	{
		public List<ParticleProcessor<S>> processors = new ArrayList<ParticleProcessor<S>>();

		@SuppressWarnings("unchecked")
		public ForkedProcessor(ParticleProcessor<S>... items)
		{
			processors = new ArrayList(Arrays.asList(items));
		}

		public ForkedProcessor(Collection<ParticleProcessor<S>> items)
		{
			this.processors = CollUtils.list(items);
		}

		public void process(S state, double weight)
		{
			for (ParticleProcessor<S> processor : processors)
				processor.process(state, weight);
		}
	}

	/**
	 * Given a function f and a loss l, returns min_x E l(f(x), f(X)), where the
	 * expectation is approximated using the samples from the SMC, and also
	 * min_x is restricted to the set of samples returned by SMC
	 * 
	 * @author bouchard
	 * 
	 * @param <D>
	 * @param <I>
	 *            I should have .equals() implemented e.g. D =
	 *            PartialCoalescentState, I = Set <Set < Languages>> i.e. clade
	 *            representation of the forest
	 */
	public static class ParticleMapperProcessor<D, I> implements
	ParticleProcessor<D> {
		private final Fct<D, I> prj;
		private final Counter<I> counter = new Counter<I>();

		public static <S> ParticleMapperProcessor<S, S> saveParticlesProcessor() {
			return new ParticleMapperProcessor<S, S>(new Id<S>());
		}

		public static ParticleMapperProcessor<PartialCoalescentState, PartialCoalescentState> saveCoalescentParticlesProcessor()
		{
			return saveParticlesProcessor();
		}

		public ParticleMapperProcessor(Fct<D, I> prj) {
			this.prj = prj;
		}

		public void process(D state, double weight) {
			counter.incrementCount(prj.evalAt(state), weight);
		}

		public I centroid(LossFct<I> loss) {
			return new BayesRiskMinimizer<I>(loss).findMin(counter);
		}

		public String printweights() {
			String s = "";
			for (I key : counter) {
				double value = counter.getCount(key);
				s += value + ",";
			}
			return s;
		}

		public I map() {
			return counter.argMax();
		}

		public I sample(Random rand) {
			double[] probs = new double[counter.size()];
			ArrayList<I> states = new ArrayList<I>();
			int i = 0;
			for (I key : counter.keySet()) {
				states.add(key);
				probs[i++] = counter.getCount(key);
			}
			int j = Multinomial.sample(rand, probs);
			return states.get(j);
		}

		public Counter<I> getCounter() {
			return counter;
		}

	}
	

	public static class PCSHash implements
	ParticleProcessor<PartialCoalescentState> {
		private Hasher hasher = new Hasher();

		public void process(PartialCoalescentState state, double weight) {
			hasher.add(weight).add(state.logLikelihood())
			.add(state.topHeight())
			.add(state.getUnlabeledArbre().deepToLispString());
		}

		public int getHash() {
			return hasher.hashCode();
		}
	}

	public static class MAPDecoder<D> implements ParticleProcessor<D> {
		private D argmax = null;
		private double max = Double.NEGATIVE_INFINITY;

		public void process(D state, double weight) {
			if (weight > max) {
				argmax = state;
				max = weight;
			}
		}

		public D map() {
			return argmax;
		}
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
		
		int nChains = 50;
		int nSamplesEachChain = 2000;
		double alpha = 1.0/3.0;
		
		//List<Taxon> leaves = MSAParser.parseMSA(alignmentInputFile).taxa();
		
		//List<Taxon> leaves = msa.taxa();
		//Taxon taxa = leaves.get(8);
		System.out.println(data.nSites());
		System.out.println(msa);
		
		String msaString = msa.toString();
		
		int Ntemp = 500;
		List<Double> originalTemp = new ArrayList<Double>();
		for(int i = 0; i < Ntemp; i++)
			originalTemp.add(Math.pow((1.0*i+1)/(Ntemp*1.0), 3));

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
			leaves.put(
			lang,
			FastDiscreteModelCalculator.observation(ctmc,
					observationSite, false));
		}
		
		double loglikelihood = UnrootedTreeState.computeLogLikelihood(urt,leaves);
		
		
		
		//double loglikelihood1 = UnrootedTreeState.computeLogLikelihoodSiteTemp(urt,leaves, sites, temp);
		//UnrootedTreeState(urt, leaves, urt.getLogPrior(), 0);
		System.out.println(data.nSites());

		System.out.println(loglikelihood);
		
		subSamplingSMCsampler<UnrootedTreeState> pc = new subSamplingSMCsampler<UnrootedTreeState>();
		pc.setNewTemperature2(originalTemp, len);
		pc.setData(data);
		Pair<List<Integer>, Pair<Integer, Double>> tempSeq = pc.tempSeqProcess(pc.getOneTemperatureSeq(498));
		Pair<List<Integer>, Pair<Integer, Double>> tempSeqDiff = pc.tempSeqProcess(pc.getOneTemperatureDiffSeq(499));
		pc.N = 100;
		
		ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;

		proposalOptions.useStochasticNearestNeighborInterchangeProposal=true;				
		LinkedList<ProposalDistribution> proposalDistributions = new LinkedList<ProposalDistribution>();
		// ParticleKernel<UnrootedTreeState> ppk
		RootedTree initrt = TreeGenerators.sampleExpNonclock(randTree, nTaxa, treeRate);
	    UnrootedTree initTree = UnrootedTree.fromRooted(rt);
		
		Gamma exponentialPrior = Gamma.exponential(10.0);
		StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
				exponentialPrior);
		UnrootedTreeState ncts = UnrootedTreeState.initFastState(initTree, data, ctmc, priorDensity);
		subSamplingKernel ppk = new subSamplingKernel(ncts, Ntemp, proposalDistributions, proposalOptions);				
		TreeDistancesProcessor tdp = new TreeDistancesProcessor();
		pc.sample(ppk, tdp);
		
		System.out.println(originalTemp.get(499));
		System.out.println(Arrays.toString(pc.getOneTemperatureSeq(499)));
		System.out.println(Arrays.toString(pc.getOneTemperatureSeq(0)));
		System.out.println(Arrays.toString(pc.getOneTemperatureDiffSeq(1)));
		//System.out.println(Arrays.toString(tempSeq.getFirst()));
		
		//modify FastDiscreteModelCalculator;





	}

}
