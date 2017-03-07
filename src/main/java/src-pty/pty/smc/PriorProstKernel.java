package pty.smc;
import java.util.*;


import nuts.math.Sampling;
import nuts.maxent.SloppyMath;

import fig.basic.Pair;
import fig.prob.Multinomial;



public class PriorProstKernel implements ParticleKernel<PartialCoalescentState>{
	private PartialCoalescentState initial;
	
	public PriorProstKernel (PartialCoalescentState initial) {
		this.initial = initial;
		if (!initial.isClock())
		  throw new RuntimeException("Not yet supported");
	}
	
	public Pair<PartialCoalescentState, Double> next (Random rand, PartialCoalescentState state) {
		final double delta = Sampling.sampleExponential (rand, 1.0/PriorPriorKernel.nChoose2(state.nRoots()));
		
		int nroots = state.nRoots();
		
		ArrayList<ArrayList<PartialCoalescentState>> candidateStates = new ArrayList<ArrayList<PartialCoalescentState>>();
		ArrayList<ArrayList<Double>> candidateProbabilities  = new ArrayList<ArrayList<Double>> ();
		
		// We need to sample a pair (i,j).
		// This is organized to avoid having to convert (i,j) to a linear order.
		// We do this by organizing each pair (i,j) by s = i+j (called firstLevel).
		// Given s, then sample i (called secondLevel)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            . 
		
		for (int i = 0 ; i < 2*(nroots-1) - 1; i++) {
			candidateStates.add (new ArrayList<PartialCoalescentState>());
			candidateProbabilities.add(new ArrayList<Double> ());
		}
		
		double[] tmparray = new double[nroots*(nroots-1)/2];
		int count = 0 ;
		for (int i = 0 ; i  < nroots; i++) {			
			for (int j = 0; j < i ; j++) {				
				PartialCoalescentState oldState = state;
				PartialCoalescentState newState = oldState.coalesce (i,j, delta, 0, 0);
				ArrayList<PartialCoalescentState> tmpState = candidateStates.get(i+j-1);
				tmpState.add (newState);
				ArrayList<Double> tmp  = candidateProbabilities.get(i+j-1);
				tmp.add(newState.logLikelihoodRatio());
				tmparray[count++] = newState.logLikelihoodRatio();
				
			}
		}
		
		double lognorm = SloppyMath.logAdd( tmparray );
		
		double[] firstLevelProbs = new double[candidateProbabilities.size()];
		for (int i = 0; i < candidateProbabilities.size(); i++) {
			ArrayList<Double> tmp = candidateProbabilities.get(i);
			firstLevelProbs[i] = Math.exp(logAdd (tmp) - lognorm);
		}
		int firstLevel = Multinomial.sample( rand, firstLevelProbs);
		ArrayList<Double> tmp =  candidateProbabilities.get(firstLevel);
		double[] secondLevelProbs = new double[tmp.size()];
		for (int i = 0 ; i < secondLevelProbs.length; i++)
			secondLevelProbs [i] = Math.exp(tmp.get(i) - lognorm)/firstLevelProbs[firstLevel];
		int secondLevel = Multinomial.sample(rand, secondLevelProbs);
		
		ArrayList<PartialCoalescentState> tmpState = candidateStates.get(firstLevel);
		PartialCoalescentState result = tmpState.get(secondLevel);
		return Pair.makePair(result, lognorm);
		
		
	}
	
	public int nIterationsLeft (PartialCoalescentState state) {
		return state.nIterationsLeft();
	}
	
	public PartialCoalescentState getInitial () {
		return initial;
	}
	
	private double logAdd (ArrayList<Double> tmp) {
		double[] tmparray = new double[tmp.size()];
		for (int i = 0 ; i < tmp.size(); i++)
			tmparray[i] = tmp.get(i);
		return SloppyMath.logAdd( tmparray) ;
		
	}
	}
