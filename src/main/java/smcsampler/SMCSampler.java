package smcsampler;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;
import monaco.process.ProcessSchedule;
import monaco.process.ProcessScheduleContext;
import monaco.process.ResampleStatus;
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
import pty.mcmc.UnrootedTreeState;
import pty.smc.PartialCoalescentState;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.Parallelizer;
import fig.prob.Multinomial;
import goblin.BayesRiskMinimizer;
import goblin.BayesRiskMinimizer.LossFct;

public final class SMCSampler<S>
{
	@Option
	public boolean verbose = true;
	public int N = 100;
	@Option
	public Random rand = new Random(1);
	@Option
	public boolean resampleLastRound = true;
	@Option
	public int nThreads = 1;
	@Option
	public ResamplingStrategy resamplingStrategy = ResamplingStrategy.ESS;

	@Option
	public  AdaptiveScheme adaptiveScheme=AdaptiveScheme.CONSTANT; 
	public static double essRatioThreshold = 0.5;

	private List<S> conditional = null; // exclude initial state
	private double[] conditionalUnnormWeights = null;

	private ProcessSchedule schedule = null;
	private double ess = N;
	private double tempDiff;

	@Option
	public boolean adaptiveTempDiff = false;
	private boolean isLastIter = false;

	public static double alpha = 0.90;

	public int adaptiveType = 0;
	public PrintWriter smcSamplerOut = null;

	public double getEss() {
		return ess;
	}

	public void setEss(double ess) {
		this.ess = ess;
	}

	public void setProcessSchedule(ProcessSchedule schedule) {
		this.schedule = schedule;
	}

	public void setConditional(List<S> conditional,
			double[] conditionalUnnormWeights)
	{
		// no I think it's ok
		// if (true)
		// throw new
		// RuntimeException("seems like conditional not used in resampling step?");
		if (conditional.size() != conditionalUnnormWeights.length)
			throw new RuntimeException();
		this.conditional = conditional;
		this.conditionalUnnormWeights = conditionalUnnormWeights;
	}

	public void setUnconditional() {
		this.conditional = null;
		this.conditionalUnnormWeights = null;
	}

	public boolean isConditional() {
		return conditional != null;
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
		for (double w : ws)
			sumOfSqr += w * w;
		return 1.0 / sumOfSqr;
	}

	public static <S> void bootstrapFilter(final SMCSamplerKernel<S> kernel,
			//			final ParticleProcessor<S> 
			final TreeDistancesProcessor processor, final int N,
			final Random rand) {
		SMCSampler<S> PF = new SMCSampler<S>();
		PF.N = N;
		PF.rand = rand;
		PF.resampleLastRound = false;
		PF.resamplingStrategy = ResamplingStrategy.ALWAYS;
		PF.sample(kernel, processor);
	}

	private long[] seeds;
	private List<S> samples;
	private double[] logWeights;
	private double[] incrementalLogWeights;
	private double varLogZ = 0;

	public List<S> getSamples() {
		return samples;
	}

	public double[] getLogWeights() {
		return logWeights;
	}

	private void propagateAndComputeWeights(final SMCSamplerKernel<S> kernel,
			final int t)
	{
		if (verbose)
			LogInfo.track("Processing...", false);

		final double[] normalizedWeights0 = logWeights.clone();
		NumUtils.expNormalize(normalizedWeights0);
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
				if (x == 0 && isConditional()) {
					samples.set(x, conditional.get(t));
					logWeights[x] = conditionalUnnormWeights[t];
				} else {
					final Pair<S, Double> current = kernel.next(rand,
							samples.get(x));
					if (current == null) {
						samples.set(x, null);
						logWeights[x] = Double.NEGATIVE_INFINITY;
					} else {
						samples.set(x, current.getFirst());
						incrementalLogWeights[x] = current.getSecond();
						logWeights[x] += incrementalLogWeights[x];
					}
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

	private void init(final SMCSamplerKernel<S> kernel)
	{
		// ancestors = new ArrayList<Integer>(N);
		// for (int n = 0; n < N; n++)
		// ancestors.add(n);
		lognorm = 0.0;
		samples = new ArrayList<S>(N);
		S initial = kernel.getInitial();
		for (int n = 0; n < N; n++)
			samples.add(initial);
		logWeights = new double[N]; // init to zero
		incrementalLogWeights = new double[N]; // init to zero
	}

	private void newProcess(int t, double[] normalizedWeights,
			//			ParticleProcessor<S> 
			TreeDistancesProcessor processor, int T) {
		if (schedule != null
				&& schedule.shouldProcess(new ProcessScheduleContext(t,
						t == T - 1, ResampleStatus.NA))) {
			if (verbose)
				LogInfo.track("Processing particles");
			for (int n = 0; n < normalizedWeights.length; n++) {
				if (verbose)
					LogInfo.logs("Particle " + (n + 1) + "/"
							+ normalizedWeights.length);
				if (samples.get(n) != null) {
					processor.process(samples.get(n), normalizedWeights[n]);
					samples.set(n, null);
				}
			}
			if (verbose)
				LogInfo.end_track();
		}
	}


	public double temperatureDifference(final double alpha,double absoluteAccuracy,double min,double max) {
		final double[] logLikePrior = new double[samples.size()];
		for (int n = 0; n < samples.size(); n++) {
			UnrootedTreeState urt = ((UnrootedTreeState) samples.get(n));
			logLikePrior[n] = urt.getLogLikelihood(); // + urt.getLogPrior();
		}
		int maxEval = 100;
		UnivariateFunction f = new UnivariateFunction() {
			public double value(double x) {
				double[] logWeightLikePriorVec = new double[samples.size()];
				for (int n = 0; n < samples.size(); n++) 
					logWeightLikePriorVec[n] = logWeights[n] + logLikePrior[n]* x;				
				NumUtils.expNormalize(logWeightLikePriorVec);
				return (ess(logWeightLikePriorVec)/logWeightLikePriorVec.length - alpha);
			}
		};
		double result = 0;
		try {
			final double relativeAccuracy = absoluteAccuracy * 0.0001;
			// final double absoluteAccuracy0 = 1.0e-6;
			PegasusSolver solver = new PegasusSolver(relativeAccuracy,
					absoluteAccuracy);
			result = solver.solve(maxEval, (UnivariateFunction) f, min, max);
		} catch (RuntimeException e) {
			LogInfo.logsForce("Solver Fail!");
			result = -1;
		}	
		return result;
	}


	/**
	 * @param <S>
	 * @param kernel
	 *            model
	 * @param tdp
	 *            what to do with the produced sample
	 */
	public void sample(final SMCSamplerKernel<S> kernel,
			final TreeDistancesProcessor tdp) {
		//			final ParticleProcessor<S> tdp){
		init(kernel);
		//        int T = kernel.nIterationsLeft(kernel.getInitial());
		//		int T = kernel.nIterationsLeft(kernel.getInitial()); //iter 0 is for initialization
		//		if (isConditional() && conditional.size() != T)
		//			throw new RuntimeException();
		if(smcSamplerOut!=null)smcSamplerOut.println(CSV.header("t", "ESS", "tempDiff"));
		double alpha0 = alpha;
		//	for (int t = 0; t < T && !isLastIter; t++) {
		int t=0;
		while(!kernel.isLastIter()){			
			kernel.setCurrentIter(t);
//			if (t > 0)  kernel.setInitializing(false);
			if (adaptiveTempDiff) {
				tempDiff = 0;
				alpha0 =adaptiveScheme.alpha(alpha0, t); 
				if (t > 0) {						
					tempDiff = temperatureDifference(
							alpha0 * ess / samples.size(), 1.0e-7, 0, 0.1);						
					if(tempDiff == -1) tempDiff = kernel.getDefaultTemperatureDifference();
				}

			} else {
				tempDiff = kernel.getDefaultTemperatureDifference();  //1.0 / (T-1);
			}
			kernel.setTemperatureDifference(tempDiff);
			if (kernel.isLastIter())
				isLastIter = true;
			if (verbose)
				LogInfo.track("Particle generation " + (t + 1) , true); //"+ "/" + T, true);
			//			if(t+1==7){ 
			//				System.out.println("Debugging");			
			//			}
			propagateAndComputeWeights(kernel, t);
			double[] normalizedWeights = logWeights.clone();
			NumUtils.expNormalize(normalizedWeights);
			if (verbose)
				LogInfo.logs("LargestNormalizedWeights="
						+ ArrayUtils.max(normalizedWeights));
			ess = ess(normalizedWeights);
			if(smcSamplerOut!=null)
			{
				smcSamplerOut.println(CSV.body(t, ess, tempDiff));
				smcSamplerOut.flush();
			}

			if (verbose)
				LogInfo.logs("RelativeESS=" + ess
						/ normalizedWeights.length);

			//			newProcess(t, normalizedWeights, tdp, T); //TODO: DOUBLE CHECK THIS!!

			if (!kernel.isLastIter()
					&& (hasNulls(samples) || resamplingStrategy
							.needResample(normalizedWeights))) {
				samples = resample(samples, normalizedWeights, rand);
				lognorm += SloppyMath.logAdd(logWeights) - Math.log(N);
				logWeights = new double[N]; // reset weights
				ess = N; // TODO: CHECK THIS!			
			}
			if (verbose)
				LogInfo.end_track();
			if ((kernel.isLastIter() || isLastIter) && schedule != null) {
				lognorm += SloppyMath.logAdd(logWeights) - Math.log(N);
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
//			if (schedule != null)
//				schedule.monitor(new ProcessScheduleContext(t, t == T - 1,
//				ResampleStatus.NA));
			t++;
		}
		setUnconditional();
		if(smcSamplerOut!=null)smcSamplerOut.close();
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

		// Counter<Integer> newSamples = new Counter<Integer>();
		// for (int n = 0; n < samples.size(); n++)
		// {
		// int index;
		// if (n == 0 && isConditional())
		// index = 0;
		// else
		// {
		// try { index = SampleUtils.sampleMultinomial(rand, weights); }
		// catch (RuntimeException re) { throw new MeasureZeroException(); }
		// }
		// newSamples.incrementCount(index, 1.0);
		// }
		// newSamples.normalize();
		// double [] resultWeight = new double[newSamples.size()];
		// List<S> resultList = new ArrayList<S>(newSamples.size());
		// int i = 0;
		// for (Integer key : newSamples.keySet())
		// {
		// resultList.add(samples.get(key));
		// resultWeight[i++] = newSamples.getCount(key);
		// }
		// return Pair.makePair(resultList,resultWeight);
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

		// List<Integer> newAncestors = new ArrayList<Integer>(N);
		// final List<S> result = new ArrayList<S>(list.size());
		// Set<Integer> indices = new HashSet<Integer>();
		// for (int n = 0; n < list.size(); n++)
		// {
		// int index;
		// if (n == 0 && isConditional())
		// index = 0;
		// else
		// {
		// try {
		// index = SampleUtils.sampleMultinomial(rand, w); }
		// catch (RuntimeException re)
		// {
		// throw new MeasureZeroException();
		// }
		// }
		// newAncestors.add(ancestors.get(index));
		// indices.add(index);
		// result.add(list.get(index));
		// }
		// if (verbose) LogInfo.logs("EffectiveSamplingSize="+indices.size());
		// ancestors = newAncestors;
		// return result;
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

}