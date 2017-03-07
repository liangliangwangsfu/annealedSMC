package pty.smc;
import java.util.List;
import java.util.Random;

import nuts.math.Sampling;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import fig.basic.Option;
import fig.basic.Pair;


public class NCPriorPriorKernel implements LazyParticleKernel<PartialCoalescentState>
{
  @Option public static double deltaProposalRate = -1.0;
  @Option public static boolean useOptimal = false;
  
  private PartialCoalescentState initial;
  
  public NCPriorPriorKernel(PartialCoalescentState initial) { setInitial(initial); }
  
  public void setInitial(PartialCoalescentState initial) 
  { 
    this.initial = initial; 
  }
  
  @Override
  public PartialCoalescentState getInitial() { return initial; }
  
  @Override
  public int nIterationsLeft(PartialCoalescentState partialState)
  {
    return partialState.nIterationsLeft();
  }
  
  @SuppressWarnings({ "unchecked", "rawtypes" })
  @Override
  public Pair<PartialCoalescentState, Double> next(Random rand,
      PartialCoalescentState current)
  {
    return (Pair) _next(rand, current, false);
  }
  public Object _next(Random rand,
      PartialCoalescentState current, boolean isPeek)
  {	  	  
    // 1-sample branch lengths
//	double rate=deltaProposalRate*PriorPriorKernel.nChoose2(current.nRoots());
	double rate=deltaProposalRate;

    final double delta0 = Sampling.sampleExponential(rand, 1.0/rate)/(current.nRoots()==2 ? 2.0 : 1.0),
                 delta1 = Sampling.sampleExponential(rand, 1.0/rate)/(current.nRoots()==2 ? 2.0 : 1.0);
    
    //System.out.println("two sampled branches: "+delta0+","+delta1);
    
     // 2- sample a random pair (without replacement)
    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, current.nRoots(), 2);
    if (sampledIndices.size() != 2)
      throw new RuntimeException();
    final int i0 = sampledIndices.get(0),
              i1 = sampledIndices.get(1);
    
    if (isPeek)
    {
      if (useOptimal) throw new RuntimeException();
      return current.peekLogLikelihoodRatio(i0, i1, 0.0, delta0, delta1) - Math.log(current.peekNNonTrivialRoots(i0, i1)); // slightly off 
//      return current.peekLogLikelihoodRatio(i0, i1, 0.0, delta0, delta1); 
    }
    else
    {
      PartialCoalescentState _next = current.coalesce(i0, i1, 0.0, delta0, delta1);
      PartialCoalescentState next = (PartialCoalescentState) _next;
      double logW = (useOptimal ? next.nonClockLogWeight() : next.logLikelihoodRatio() - Math.log(next.nNonTrivialRoots()));
//      double logW = (useOptimal ? next.nonClockLogWeight() : next.logLikelihoodRatio());
      return Pair.makePair(_next, logW);
    }
  }

  @Override
  public double peekNext(Random rand, PartialCoalescentState current)
  {
    return (Double) _next(rand, current, true);
  }
}