package pty.smc;
import java.io.*;
import java.util.*;

import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.ParticleFilter.ResamplingStrategy;

import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.Parallelizer;
import monaco.process.ProcessSchedule;
import nuts.lang.ArrayUtils;
import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils.*;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.MathUtils;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public final class LazyParticleFilter<T>
{
	private final LazyParticleKernel<T> kernel;
	private final ParticleFilterOptions options;
	//  private final ProcessSchedule schedule; 

	public LazyParticleFilter(LazyParticleKernel<T> kernel,
			ParticleFilterOptions options)
	{
		this.kernel = kernel;
		this.options = options;
	}

	public static class ParticleFilterOptions
	{
		@Option public boolean verbose = false;
		@Option public int nParticles = 100000;
		@Option public int maxNUniqueParticles = 1000; // memory bottleneck
		@Option public Random rand = new Random(1);
		@Option public boolean resampleLastRound = true;
		@Option public int nThreads = 8;
		@Option public int processBatchSize = 1000;
		@Option public int finalMaxNUniqueParticles = Integer.MAX_VALUE;
		@Option public boolean parallelizeFinalParticleProcessing = false;
		@Option public double populationShrinkFactor = 0.9;

		public void check()
		{
			if (populationShrinkFactor >= 1)
				throw new RuntimeException();
			if (maxNUniqueParticles <= 1)
				throw new RuntimeException();
		}
	}

	public static class Eager2LazyAdaptor<S> implements LazyParticleKernel<S>
	{
		private final ParticleKernel<S> kernel;
		public Eager2LazyAdaptor(ParticleKernel<S> kernel)
		{
			this.kernel = kernel;
		}
		@Override
		public Pair<S, Double> next(Random rand, S current)
		{
			return kernel.next(rand, current);
		}
		@Override
		public int nIterationsLeft(S partialState)
		{
			return kernel.nIterationsLeft(partialState);
		}
		@Override
		public S getInitial()
		{
			return kernel.getInitial();
		}
		@Override
		public double peekNext(Random rand, S current)
		{
			return kernel.next(rand, current).getSecond();
		}
	}

	public static interface LazyParticleKernel<S> extends ParticleKernel<S>
	{	
		public double peekNext(Random rand, S current);
	}

	public static interface Distribution<T>
	{
		public List<T> sampleNTimes(Random r, int n);
	}

	private ParticlePopulation pruneIfNeeded(Random rand, ParticlePopulation result, int maxNUniqueParticles)
	{
		double currentExpectedSupportSize = result.expectedResampledSupportSize(options.nParticles);
		if (currentExpectedSupportSize > maxNUniqueParticles)
		{
			ParticlePopulation resampled = null;
			double nResampled = result.nParticlesTokens();
			while (currentExpectedSupportSize > maxNUniqueParticles)
			{
				resampled = result.prunePopulation(rand, 1+(int) nResampled);
				currentExpectedSupportSize = resampled.expectedResampledSupportSize(options.nParticles);
				nResampled = options.populationShrinkFactor  * nResampled;
			}
			logs("Current population pruned to " + (1+(int) nResampled) + " (" + resampled.nParticlesTokens() + " unique particles)");
			result = resampled;
		}
		return result;
	}

	public Pair<ParticlePopulation,Double> extend(Random rand, final List<T> parentParticles)
	{
		final int size = parentParticles.size();
		// create seeds
		final long [] seeds = new long[size];
		for (int i = 0; i < seeds.length; i++)
			seeds[i] = rand.nextLong();

		// container to handle result
		final double [] logWeights = new double[size];

		track("Extension..");
		// do the jobs
		Parallelizer<Integer> parallelizer = new Parallelizer<Integer>(options.nThreads);
		parallelizer.setPrimaryThread();
		parallelizer.process(CollUtils.ints(size), new Parallelizer.Processor<Integer>() {
			public void process(Integer x, int _i, int _n, boolean log) {
				final T parent = parentParticles.get(x);
				final Random curRand = new Random(seeds[x]);
				final double currentW = kernel.peekNext(curRand, parent);
				logWeights[x] = currentW;
			}});
		end_track();

		// log normalizer
		double logSum = Double.NEGATIVE_INFINITY;
		for (final double n : logWeights)
			logSum = SloppyMath.logAdd(logSum, n);
		//    for(int i=0;i<logWeights.length;i++) System.out.print(logWeights[i]+" ");
		//   System.out.println();
		NumUtils.expNormalize(logWeights);
		Object [] currents = new Object[logWeights.length];
		ParticlePopulation newPop = new ParticlePopulation(logWeights, parentParticles.toArray(), seeds, currents);
		return Pair.makePair(newPop, logSum);
	}



	public double sample(ParticleProcessor<T> ...particleProcessors)
	{
		track("Sampling PF");
		Random curRandom = null;
		synchronized (this) { curRandom = new Random(options.rand.nextLong()); }
		Distribution<T> previousPopulation = initialPopulation;
		double zHat = 0.0;
		final int T = kernel.nIterationsLeft(kernel.getInitial());
		ParticlePopulation finalPopulation = null;
		for (int t = 0; t < T; t++)
		{
			track("Generation " + t + "/" + T);

			Pair<ParticlePopulation,Double> resultPair = extend(curRandom, previousPopulation.sampleNTimes(curRandom, options.nParticles));
			ParticlePopulation newPopulation = resultPair.getFirst();
			logs("maxNormalizedWeight=" + newPopulation.maxNormalizedWeight());
			newPopulation = pruneIfNeeded(curRandom, newPopulation, options.maxNUniqueParticles);
			logs("expectedResampledSupportSize=" + newPopulation.expectedResampledSupportSize(options.nParticles));

			zHat += resultPair.getSecond() - Math.log(options.nParticles);
			previousPopulation = newPopulation;    
			if (t == T-1)
				finalPopulation = newPopulation;
			end_track();
		}
		end_track();

		if (options.resampleLastRound)
			finalPopulation = pruneIfNeeded(curRandom, finalPopulation, options.finalMaxNUniqueParticles );

		finalPopulation.process(particleProcessors);

		return zHat;
	}


	private final Distribution<T> initialPopulation = new Distribution<T>() {
		@Override
		public List<T> sampleNTimes(Random r, int n)
		{
			final T initial = kernel.getInitial();
			List<T> result = list();
			for (int i = 0; i < n; i++)
				result.add(initial);
			return result;
		}
	};

	public class  ParticlePopulation implements Distribution<T>
	{
		private final double [] normalizedWeights;
		private final Object [] parents, currents;
		private final long [] randomSeeds;

		private ParticlePopulation(
				double [] normalizedWeights,
				Object [] parents,
				long [] randomSeeds,
				Object [] currents)
		{
			final int size = normalizedWeights.length;
			if (parents.length != size || randomSeeds.length != size || currents.length != size)
				throw new RuntimeException();
			this.normalizedWeights = normalizedWeights;
			MathUtils.checkIsProb(normalizedWeights);
			this.parents = parents;
			this.currents = currents;
			this.randomSeeds = randomSeeds;
		}

		public double maxNormalizedWeight()
		{
			return ArrayUtils.max(normalizedWeights);
		}

		public ParticlePopulation prunePopulation(Random r, int nResampling)
		{
			Counter<Integer> indices = Sampling.efficientMultinomialSampling(r, normalizedWeights, nResampling);
			double [] normalizedWeights = new double[indices.keySet().size()];
			Object [] parents = new Object[indices.keySet().size()];
			long [] randomSeeds = new long[indices.keySet().size()];
			Object [] currents = new Object[indices.keySet().size()];
			int i = 0;
			for (int index : indices)
			{
				normalizedWeights[i] = indices.getCount(index);
				parents[i] = this.parents[index];
				randomSeeds[i] = this.randomSeeds[index];
				currents[i] = this.currents[index];
				if (parents[i] == null && currents[i] == null)
					throw new RuntimeException();
				i++;
			}
			NumUtils.normalize(normalizedWeights);
			return new ParticlePopulation(normalizedWeights, parents, randomSeeds, currents);
		}

		public void process(final ParticleProcessor<T>[] particleProcessors)
		{
			if (particleProcessors.length == 0)
				return;

			Set<Integer> allIndices = set(CollUtils.ints(nParticlesTokens()));
			track("Processing particles");
			while (!allIndices.isEmpty())
			{
				// divide by blocks to avoid to store too much in memory!
				final Set<Integer> subset = set();
				int i = 0;
				forLoop:for (int index : allIndices)
				{
					subset.add(index);
					if (i++ >  options.processBatchSize)
						break forLoop;
				}
				allIndices.removeAll(subset);
				unLazyParticles(subset);
				int nThreads = options.parallelizeFinalParticleProcessing  ? options.nThreads : 1;
				Parallelizer<Integer> parallelizer = new Parallelizer<Integer>(nThreads);
				parallelizer.setPrimaryThread();
				parallelizer.process(list(subset), new Parallelizer.Processor<Integer>() {
					public void process(Integer x, int _i, int _n, boolean log) {
						final T state = (T) currents[x];
						final double w = normalizedWeights[x];
						if (log)
							logs("" + (_i) + "/" + _n);
						for (ParticleProcessor<T> pro : particleProcessors)
							pro.process(state, w);
						// de-assign to save memory
						currents[x] = null;
					}});
			}
			end_track();
		}

		private void unLazyParticles(Set<Integer> indices)
		{
			track("Unlazy...");
			List<Integer> listOfUniqueIndices = list();
			for (Integer index : indices)
				if (currents[index] == null)
					listOfUniqueIndices.add(index);
			Parallelizer<Integer> parallelizer = new Parallelizer<Integer>(options.nThreads);
			parallelizer.setPrimaryThread();
			parallelizer.process(listOfUniqueIndices, new Parallelizer.Processor<Integer>() {
				public void process(Integer x, int _i, int _n, boolean log) {
					final Random rand = new Random(randomSeeds[x]);
					final T parent = (T) parents[x];
					final T current = kernel.next(rand, parent).getFirst();
					parents[x] = null;
					randomSeeds[x] = -1;
					currents[x] = current;
				}});
			end_track();
		}

		public int nParticlesTokens()
		{
			return normalizedWeights.length;
		}


		public double expectedResampledSupportSize(int nResampling)
		{
			return Sampling.expectedNumberDistinctParticles(normalizedWeights, nResampling);
		}

		public List<T> sampleNTimes(Random r, int n)
		{
			// sample indices
			Counter<Integer> indices = Sampling.efficientMultinomialSampling(r, normalizedWeights, n);
			// un-lazy the entries that have not been un-lazied already
			unLazyParticles(indices.keySet());      
			List<T> result = list();
			for (Integer index : indices.keySet())
			{
				T current = (T)currents[index];
				for (int i = 0; i < indices.getCount(index); i++)
					result.add(current);
			}      
			return result;
		}

	}
	private void logs(Object s)
	{
		if (options.verbose)
			LogInfo.logs(s);
	}
	private void track(Object s)
	{
		if (options.verbose)
			LogInfo.track(s);
	}
	private void end_track()
	{
		if (options.verbose)
			LogInfo.end_track();
	}

}
