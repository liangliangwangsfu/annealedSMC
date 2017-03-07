package pty.smc;
import java.util.*;

import nuts.math.Sampling;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;

import pty.smc.LazyParticleFilter.LazyParticleKernel;

/**
 * See Yee Whye Teh et al. Bayesian Agglomerative Clustering with Coalescents
 * @author bouchard
 */
public class PriorPriorKernel implements LazyParticleKernel<PartialCoalescentState>
{
  @Option public static boolean printBranchLengthMagnitudes = false;
  private final PartialCoalescentState initial;
  public PriorPriorKernel(PartialCoalescentState initial) { this.initial = initial; }
  public PartialCoalescentState getInitial() { return initial; }
  public int nIterationsLeft(PartialCoalescentState partialState)
  {
    return partialState.nIterationsLeft();
  }
  
//  public PartialCoalescentState coalesce(int left, int right, 
//      double delta, double leftIncrement, double rightIncrement) 
  
  
  public Object _next(
      Random rand,
      PartialCoalescentState current, boolean isPeek)
  {	
    // 1- sample the exponential waiting time
    //final double delta = Sampling.sampleExponential(rand, (initial.isClock() ? 1.0 : 0.5)/nChoose2(current.nRoots()));        
	  final double delta = Sampling.sampleExponential(rand, (initial.isClock() ? 0.1 : 0.05)/nChoose2(current.nRoots()));   
    
    // 2- sample a random pair (without replacement)
    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, current.nRoots(), 2);
    	        
    final int i0 = sampledIndices.get(0),
              i1 = sampledIndices.get(1);    
    		    		
    double leftIncrement = 0.0,
           rightIncrement= 0.0;
    if (!initial.isClock())
    {
      double incr = Sampling.sampleExponential(rand, 0.1);
      if (rand.nextBoolean())
        leftIncrement += incr;
      else
        rightIncrement+= incr;
    }
    
    
    if (printBranchLengthMagnitudes)
    {
      LogInfo.logsForce("Merged things:" + current.clade(i0)+ "," + current.clade(i1));
      
      // show likelihood for different times
      for (int i = -20; i < 20; i += 2)
      {
        PartialCoalescentState result = current.coalesce(
            i0, i1, 
            delta * Math.pow(2.0, i), leftIncrement, rightIncrement);
        LogInfo.logsForce("\t"+ i + "\t"+ delta * Math.pow(2.0, i) + "\t" + result.logLikelihoodRatio());
      }
    }

//    if (returnProposalLogRatio)
//    {
//      PhyloParticle result = null;
//      if (!isPeek)
//        result = current.coalesce(
//            i0, i1, 
//            delta, leftIncrement, rightIncrement);
//      return Pair.makePair(result, second)
//    }
        

    PartialCoalescentState result = null;
    Double logw = null;
    if (isPeek)
      logw = current.peekLogLikelihoodRatio(i0, i1, delta, leftIncrement, rightIncrement);
    else
      result = current.coalesce(
        i0, i1, 
        delta, leftIncrement, rightIncrement);
    
    // 
//    System.out.println("Current Log likelihood=" + result.logLikelihood());
//    System.out.println("Current Log likelihood ratio=" + result.logLikelihoodRatio());
//    
//    {  // TODO: REMOVE ME!!!
//      LogInfo.logsForce("Current delta:" + delta);
//      for (int i = 0; i < current.nRoots(); i++)
//      {
//        Set<Language> s1 = current.clade(i);
//        LogInfo.logsForce("Cost for " + s1 + " and...");
//        Counter<Integer> dist = new Counter<Integer>(),
//         fulldist = new Counter<Integer>();
//        for (int j = 0; j < current.nRoots(); j++)
//          if (i != j)
//          {
//            
//            Set<Language> s2 = current.clade(j);
//            PartialCoalescentState temp = current.coalesce(
//                i,j, 
//                delta+current.topHeight());
//            dist.setCount(j, temp.logLikelihoodRatio());
//            fulldist.setCount(j, temp.logLikelihood());
//          }
//        for (int key : dist)
//        {
//          Set<Language> otherClade = current.clade(key);
//          LogInfo.logsForce("\t" + otherClade + ":" + dist.getCount(key) + " (full="+ fulldist.getCount(key) + ")");
//        }
//      }
//    }
    
//    {// TODO: remove me
//      Set<Language> s1 = current.clade(sampledIndices.get(0));
//      Set<Language> s2 = current.clade(sampledIndices.get(1));
//      LogInfo.logsForce("Initial:" + SymmetricDiff.clades2arbre(current.allClades()).deepToString());
//      LogInfo.logsForce("Coalescing clades " + sampledIndices.get(0) 
//          + "(" + s1.toString ()+ ")," + sampledIndices.get(1)  + "(" + s2.toString() + ")" + "\tweight = "  + result.lastStepLogLikelihoodRatio(), true);
//      LogInfo.logsForce("Obtained:" + SymmetricDiff.clades2arbre(result.allClades()).deepToString());
//    }
    
//    {
//      /// TODO: remove me
//      if (result.isFinalState())
//      {
//        CTMC ctmc = result.getCTMC();
//        double logLikelihood = 0.0;
//        for (int s = 0; s < ctmc. nSites(); s++)
//        {
//          GMFct<Language> pots = DiscreteBP.toGraphicalModel(result, s);
////          LogInfo.logsForce("Pots:" + pots);
//          TreeSumProd<Language> tsp =new TreeSumProd<Language>(pots);
////          LogInfo.logsForce("Post:" + tsp.moments());
//          logLikelihood += tsp.logZ();
//        }
//        MathUtils.checkClose(logLikelihood, result.logLikelihood());
//      }
//    }
    
    // 3- the weight update is simply equal to the ratio of the new likelihood score to the old one
    if (isPeek)
      return logw;
    else
      return Pair.makePair(result, result.logLikelihoodRatio());
  }
  public static double nChoose2(double n) { return n*(n-1)/2; }
  
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

}
