
package pty.smc;
import java.io.*;
import java.util.*;

import pty.RootedTree;
import pty.smc.PartialCoalescentState.CoalescentNode;

import fig.basic.Option;
import fig.basic.Pair;
import goblin.Taxon;
import nuts.math.Sampling;
import nuts.util.CollUtils.*;
import nuts.util.Arbre;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;

/*
public class NCPriorPriorKernel2 implements ParticleKernel
{
//  @Option public static double deltaProposalRate = -1.0;
	@Option public static double deltaProposalRate = 1.0;
  @Option public static boolean useOptimal = true;
  @Option public static boolean useLazy = false;
  
  private Object initial;
  
  public NCPriorPriorKernel2(Object initial) { setInitial(initial); }
  
  public void setInitial(Object initial) 
  { 
    if (useLazy)
      this.initial = new LazyPCS((PartialCoalescentState) initial);
    else
      this.initial = initial; 
  }
  
  @Override
  public Object getInitial() { return initial; }
  
  @Override
  public int nIterationsLeft(Object partialState)
  {
    if (useLazy)
      return ((LazyPCS) partialState).nIterationsLeft();
    else
      return ((PartialCoalescentState) partialState).nIterationsLeft();
  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  @Override
  public Pair<PartialCoalescentState, Double> next(Random rand,
      Object _current)
  {
    PartialCoalescentState current = (useLazy ? ((LazyPCS) _current).getState() : (PartialCoalescentState) _current);
    
    // 1-sample branch lengths
    final double delta = Sampling.sampleExponential(rand, 1.0/deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0); 
    
     
//    final double delta0 = Sampling.sampleExponential(rand, 1.0/deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0),
//                 delta1 = Sampling.sampleExponential(rand, 1.0/deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0);
       
     // 2- sample a random pair (without replacement)
    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, current.nRoots(), 2);
    if (sampledIndices.size() != 2)
      throw new RuntimeException();
    final int i0 = sampledIndices.get(0),
              i1 = sampledIndices.get(1);
        
//    double currentHeight=current.maxSubHeight();
    
    double bl0=RootedTree.Util.averageHeight(current.getSubtree(i0)),bl1=RootedTree.Util.averageHeight(current.getSubtree(i1));
    double currentHeight=current.maxSubHeight();
    double b0=currentHeight+delta-bl0, b1=currentHeight+delta-bl1;
    double delta0=(b0+b1)*rand.nextDouble(), delta1=b0+b1-delta0; 
    
// System.out.println("averageHeight "+averageHeight + "delta "+ delta+ " delta0: "+ delta0+ " delta1: "+ delta1 + "current.getHeight(i0): "+RootedTree.Util.averageHeight(current.getSubtree(i0)));
    
    final Object _next = (useLazy ? new LazyPCS(current, 0.0, i0, i1, delta0, delta1) : current.coalesce(i0, i1, 0.0, delta0, delta1));
    
    // 3 - compute LOG weight..
    double logW = Double.NaN;
    if (!useLazy) 
    {
      PartialCoalescentState next = (PartialCoalescentState) _next;
      logW = (useOptimal ? next.nonClockLogWeight() : next.logLikelihoodRatio())  - Math.log(next.nNonTrivialRoots());
//      System.out.println(next.nonClockLogWeight() +", "+ next.logLikelihoodRatio());
    }
    else
    {
      if (useOptimal) throw new RuntimeException();
      logW = current.peekLogLikelihoodRatio(i0, i1, 0.0, delta0, delta1) - Math.log(current.peekNNonTrivialRoots(i0, i1)); // slightly off
    }
       logW=logW-Sampling.exponentialLogDensity(1.0/deltaProposalRate, delta)+Math.log(b0+b1)+Sampling.exponentialLogDensity(1.0/deltaProposalRate, delta0)+Sampling.exponentialLogDensity(1.0/deltaProposalRate, delta1); 
    Pair result =  Pair.makePair(_next, logW);
    return result;
  }
  
  
 

}
*/