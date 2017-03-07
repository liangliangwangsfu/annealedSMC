package pty.smc;

import java.util.*;


import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.Arbre;
import nuts.util.Counter;

import fig.basic.LogInfo;
import fig.basic.Pair;
import fig.prob.Multinomial;

import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.test.TestBrownianModel;


public class PostPostKernel implements ParticleKernel<PartialCoalescentState>{
	private PartialCoalescentState initial;
	
	
	public PostPostKernel (PartialCoalescentState initial) {
		this.initial = initial;
		if (!initial.isClock())
		  throw new RuntimeException("Not yet supported");
	}
	
	
	public Pair<PartialCoalescentState, Double> next (Random rand, PartialCoalescentState state) {
		
		int nroots = state.nRoots();
		
		ArrayList<ArrayList<PartialCoalescentState>> candidateStates = new ArrayList<ArrayList<PartialCoalescentState>>();
		ArrayList<ArrayList<Double>> candidateProbabilities  = new ArrayList<ArrayList<Double>> ();
		// We need to sample a pair (i,j).
		// This is organized to avoid having to convert (i,j) to a linear order.
		// We do this by organizing each pair (i,j) by s = i+j (called firstLevel).
		// Given s, then sample i (called secondLevel)                    
		
		for (int i = 0 ; i < 2*(nroots-1) - 1; i++) {
			candidateStates.add (new ArrayList<PartialCoalescentState>());
			candidateProbabilities.add(new ArrayList<Double> ());
		}
		
//		LogInfo.track("prs",true);
		Counter<String> mergers = new Counter();
		
		double[] tmparray = new double[nroots*(nroots-1)/2];
		int count = 0 ;
		for (int i = 0 ; i  < nroots; i++) {			
			for (int j = 0; j < i ; j++) {				
				PartialCoalescentState oldState = state;
				List<Arbre<CoalescentNode>> roots = oldState.getRoots();
				CoalescentNode  node1 = roots.get(i).getContents(),
                  				node2 = roots.get(j).getContents();
			    BrownianModelCalculator bmc = (BrownianModelCalculator)node1.likelihoodModelCache;
			    double weight = bmc.evaluatePair(node1,node2,oldState.nRoots(), oldState.topHeight());			    
				ArrayList<Double> tmp  = candidateProbabilities.get(i+j-1);
				tmp.add(weight);
				tmparray[count++] = weight;
				
				mergers.setCount("Weight for " + node1.nodeIdentifier + "," + node2.nodeIdentifier + "=" , weight);
				
			}
		}
		
		
//		for (String entry : mergers)
//		  LogInfo.logs(entry + mergers.getCount(entry));
//		LogInfo.end_track();
		
		final double lognorm = SloppyMath.logAdd( tmparray );
		double particleweight = lognorm; 
		
		double[] firstLevelProbs = new double[candidateProbabilities.size()];
		for (int i = 0; i < candidateProbabilities.size(); i++) {
			ArrayList<Double> tmp = candidateProbabilities.get(i);
			firstLevelProbs[i] = Math.exp(logAdd (tmp) - lognorm);
		}
		int firstLevel = Multinomial.sample( rand, firstLevelProbs);
		ArrayList<Double> tmp =  candidateProbabilities.get(firstLevel);
		
//		System.out.println ("Second level");
//		for (int i = 0; i < tmp.size(); i++) {
//			System.out.print(tmp.get(i) + "\t");
//		}
//		System.out.print("\n");
		
		double[] secondLevelProbs = new double[tmp.size()];
		for (int i = 0 ; i < secondLevelProbs.length; i++)
			secondLevelProbs [i] = Math.exp(tmp.get(i) - lognorm)/firstLevelProbs[firstLevel];
		int secondLevel = Multinomial.sample(rand, secondLevelProbs);
		double weight = tmp.get(secondLevel);
		
		Pair<Integer, Integer> p = level2index (firstLevel, secondLevel);
		int i = p.getFirst();
		int j = p.getSecond();
		List<Arbre<CoalescentNode>> roots = state.getRoots();
		CoalescentNode  node1 = roots.get(i).getContents(),
          				node2 = roots.get(j).getContents();
	    BrownianModelCalculator bmc = (BrownianModelCalculator)node1.likelihoodModelCache;
	    
	    Pair<Double, Double> result  = bmc.sampleBranchLength (node1,node2, state.topHeight(), state.nRoots(), weight, rand);
	    double delta = result.getFirst();
	    double correct = result.getSecond();
	    particleweight += correct;
	    
	    PartialCoalescentState newState;
	    if (TestBrownianModel.test) {	    
	    	newState = state.coalesce (0,1,1, 0, 0);
//	    	System.out.println ("Particle weight = " + particleweight);
	    }
	    
	    newState = state.coalesce (i,j, delta, 0, 0);
//		System.out.println ("Coalescing  " + node1.nodeIdentifier + "," + node2.nodeIdentifier);
	    
		return Pair.makePair(newState, particleweight);
		
		
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
	
	public Pair<Integer, Integer> level2index (int firstLevel, int secondLevel) {
		firstLevel++;
		int i,j;
		if (firstLevel%2==1) {
			j=firstLevel/2;
			i=firstLevel-j;
		} else {
			j=firstLevel/2-1;
			i=firstLevel/2+1;
		}
		
		i += secondLevel;
		j -= secondLevel;
		return new Pair<Integer, Integer> (i,j);
	}
}
