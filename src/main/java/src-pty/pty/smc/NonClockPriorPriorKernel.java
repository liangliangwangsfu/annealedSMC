package pty.smc;
import java.io.*;
import java.util.*;

import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.PartialCoalescentState.CoalescentNode;

import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import nuts.math.Sampling;
import nuts.util.CollUtils.*;
import nuts.util.Arbre;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;

/*
public class NonClockPriorPriorKernel implements LazyParticleKernel<PartialCoalescentState>
{
  @Option public static double deltaProposalRate = 1.0;
  @Option public static boolean useOptimal = true;
  @Option public static boolean useLazy = false;
  
  @Option public static boolean printBranchLengthMagnitudes = false;
  private final PartialCoalescentState initial;
  public NonClockPriorPriorKernel(PartialCoalescentState initial) { this.initial = initial; }
  public PartialCoalescentState getInitial() { return initial; }
  public int nIterationsLeft(PartialCoalescentState partialState)
  {
    return partialState.nIterationsLeft();
  }
   
  public Object _next(
	      Random rand,
	      PartialCoalescentState current, boolean isPeek)
  {    
    // 1-sample branch lengths
    final double delta0 = Sampling.sampleExponential(rand, 1.0/deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0),
                 delta1 = Sampling.sampleExponential(rand, 1.0/deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0);
    //System.out.println("branches: "+ delta0+", "+delta1);
       
     // 2- sample a random pair (without replacement)
    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, current.nRoots(), 2);
    if (sampledIndices.size() != 2)
      throw new RuntimeException();
    final int i0 = sampledIndices.get(0),
              i1 = sampledIndices.get(1);
    
    PartialCoalescentState result = null;
    Double logw = null;

    // 3- the weight update is simply equal to the ratio of the new likelihood score to the old one
    if (isPeek)
    {    	
    	logw = useOptimal ?current.peekOptimalLogLikelihoodRatio(i0, i1, 0.0, delta0, delta1):  current.peekLogLikelihoodRatio(i0, i1, 0.0, delta0, delta1)+logCountRatio(current, i0, i1);    	
    	return logw; 
    }
    else
    {
      result = current.coalesce(i0, i1, 0.0, delta0, delta1);
      logw = (useOptimal ? result.nonClockLogWeight() : result.logLikelihoodRatio() +logCountRatio(current, i0, i1));
      return Pair.makePair(result, logw);
    }     
  }

  
  
@Override
public Pair<PartialCoalescentState, Double> next(Random rand,
    PartialCoalescentState current)
{
  return (Pair) _next(rand, current, false);
}
@Override
public double peekNext(Random rand, PartialCoalescentState current)
{
  return (Double) _next(rand, current, true);
}
public static double logCountRatio(PartialCoalescentState current, int i0, int i1)
{
	  int nRoots=current.nRoots();
	  List<Arbre<CoalescentNode>> roots= current.getRoots();
	  double denom=0;
	  for(int i=0;i<nRoots;i++)
	  {
//		  System.out.println(roots.get(i).nLeaves()-1);
		  denom+=roots.get(i).nLeaves()-1;
	  }
	  double numerator=roots.get(i0).nLeaves()+roots.get(i1).nLeaves()-1;
	  if(denom==0)return 0.0;
	  return Math.log(numerator)-Math.log(denom+1);
}
public static double nChoose2(double n) { return n*(n-1)/2; }

}
*/