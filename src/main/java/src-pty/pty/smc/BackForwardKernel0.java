package pty.smc;
import java.util.List;
import java.util.Random;

import nuts.math.Sampling;
import fig.basic.Option;
import fig.basic.Pair;

/**
 * 
 * @author Liangliang Wang
 */
public class BackForwardKernel0 implements
		ParticleKernel<PartialCoalescentState4BackForwardKernel> {
	@Option
	public static boolean printBranchLengthMagnitudes = false;
	private final PartialCoalescentState4BackForwardKernel initial;

	public BackForwardKernel0(PartialCoalescentState4BackForwardKernel initial) {
		this.initial = initial;
	}

	public PartialCoalescentState4BackForwardKernel getInitial() {
		return initial;
	}

	public int nIterationsLeft(
			PartialCoalescentState4BackForwardKernel partialState) {
		return partialState.nIterationsLeft();
	}

	public Object _next0(Random rand,
			PartialCoalescentState4BackForwardKernel current,
			boolean isPeek) {

		final double delta = Sampling.sampleExponential(rand,
				(initial.isClock() ? 0.1 : 0.05) / nChoose2(current.nRoots()));

		// 2- sample a random pair (without replacement)
		List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand,
				current.nRoots(), 2);

		final int i0 = sampledIndices.get(0), i1 = sampledIndices.get(1);

		double leftIncrement = 0.0, rightIncrement = 0.0;
		if (!initial.isClock()) {
			double incr = Sampling.sampleExponential(rand, 0.1);
			if (rand.nextBoolean())
				leftIncrement += incr;
			else
				rightIncrement += incr;
		}

		PartialCoalescentState4BackForwardKernel result = null;
		Double logw = null;
		if (isPeek)
			logw = current.peekLogLikelihoodRatio(i0, i1, delta, leftIncrement,
					rightIncrement);
		else
 {
			PartialCoalescentState resultpcs = current
					.coalesce(i0, i1, delta, leftIncrement,
					rightIncrement);
			result = new PartialCoalescentState4BackForwardKernel(resultpcs,
					null, delta);
			}

		if (isPeek)
			return logw;
		else
			return Pair.makePair(result, result.logLikelihoodRatio());
	}

	public Object _next1(Random rand,
			PartialCoalescentState4BackForwardKernel current,
			boolean isPeek) {
		// 1- go back to current state's parent partial state
		PartialCoalescentState4BackForwardKernel parent = current.parentState();
		double oldLogLikelihoodRatio = current.logLikelihoodRatio();
		double deltaOld = current.getDeltaOld();
		
		// 2- sample the exponential waiting time
		
		double param0=(initial.isClock() ? 0.1 : 0.05)
		/ nChoose2(current.nRoots() + 1);
		final double delta0 = Sampling.sampleExponential(rand, param0);

		double param1 = (initial.isClock() ? 0.1 : 0.05)
				/ nChoose2(current.nRoots());
		final double delta1 = Sampling.sampleExponential(rand, param1);

		// 3- sample a random pair (without replacement)
		List<Integer> sampledIndices0 = Sampling.sampleWithoutReplacement(rand,
				current.nRoots() + 1, 2);
		
		// double proposalLogDensity = Sampling.exponentialLogDensity(param0,
		// delta0) + Sampling.exponentialLogDensity(param1, delta1);

		double logExpDensityDeltaOld = Sampling.exponentialLogDensity(param0,
				deltaOld);

		final int i00 = sampledIndices0.get(0), i01 = sampledIndices0.get(1);
		final int i10 = sampledIndices0.get(0), i11 = sampledIndices0.get(1);
		double leftIncrement0 = 0.0, rightIncrement0 = 0.0;
		double leftIncrement1 = 0.0, rightIncrement1 = 0.0;
		if (!initial.isClock()) {
			double incr0 = Sampling.sampleExponential(rand, 0.1), incr1 = Sampling
					.sampleExponential(rand, 0.1);
			if (rand.nextBoolean()) {
				leftIncrement0 += incr0;
				leftIncrement1 += incr1;
			}
 else {
				rightIncrement0 += incr0;
				rightIncrement1 += incr1;
			}
		}


		PartialCoalescentState4BackForwardKernel result0 = (PartialCoalescentState4BackForwardKernel) parent
				.coalesce(i00, i01, delta0, leftIncrement0, rightIncrement0);

		PartialCoalescentState4BackForwardKernel result = null;
		 Double loglikeRatio0=result0.logLikelihoodRatio();
		Double logw =null;
		if (isPeek) {
			logw= loglikeRatio0+result0.peekLogLikelihoodRatio(i10, i11, delta1,
							leftIncrement1, rightIncrement1)
					- oldLogLikelihoodRatio - logExpDensityDeltaOld;
		} else {
			result = (PartialCoalescentState4BackForwardKernel) result0
					.coalesce(i10, i11, delta1, leftIncrement1, rightIncrement1);
			result.setParent(result0);
			result.setDeltaOld(delta1);
		}

		
		

		// 4- the weight update is simply equal to the ratio of the new
		// likelihood score to the old one
		if (isPeek)
			return logw;
		else
			return Pair.makePair(result, result.logLikelihoodRatio()
					+ loglikeRatio0 - oldLogLikelihoodRatio
					- logExpDensityDeltaOld);
	}

	public static double nChoose2(double n) {
		return n * (n - 1) / 2;
	}


	public Object _next(Random rand,
			PartialCoalescentState4BackForwardKernel current, boolean isPeek) {
		if (current.hasParent())
			return _next1(rand, current, isPeek);
		else
			return _next0(rand, current, isPeek);
	}

	@Override
	public Pair<PartialCoalescentState4BackForwardKernel, Double> next(
			Random rand, PartialCoalescentState4BackForwardKernel current) {
		return (Pair) _next(rand, current, false);

	}



}
