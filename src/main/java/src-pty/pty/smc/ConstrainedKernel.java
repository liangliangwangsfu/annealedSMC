package pty.smc;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils;
import nuts.util.CoordinatesPacker;
import pty.eval.SymmetricDiff;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.prob.SampleUtils;
import goblin.Taxon;

public class ConstrainedKernel implements ParticleKernel<PartialCoalescentState>
{
  @Option public static boolean usesPriorPost = false;
  @Option public static double constraintViolationLogPenalty = -100;
  private final PartialCoalescentState initial;
  private final Set<Set<Taxon>> constraints;
  public ConstrainedKernel(
      PartialCoalescentState initial, 
      Set<Set<Taxon>> constraints) { 
    this.constraints = constraints;
    this.initial = initial;
    if (!initial.isClock())
      throw new RuntimeException("Not yet supported");
//    checkCompatible(initial,constraints);
//    if (!initial.isReversible())
//      throw new RuntimeException();
  }
  public static void checkCompatible(PartialCoalescentState initial,
      Set<Set<Taxon>> constraints)
  {
    Set<Taxon> 
      s1 = initial.getObservations().observations().keySet(),
      s2 = SymmetricDiff.allLeaves(constraints);
    if (!s1.equals(s2))
      LogInfo.warning("Leaves of constraints:" +
         s1 + "\nLeaves of PCS:" + s2 +
         "\nSymmetric diff:" + CollUtils.symmetricDifference(s1,s2));
  }
  public Pair<PartialCoalescentState, Double> next (Random rand, PartialCoalescentState state) {
    int nRoots = state.nRoots();
    final double delta = Sampling.sampleExponential (rand, 1.0/PriorPriorKernel.nChoose2(nRoots));
    double [] prs = new double[nRoots*nRoots];
    double lognorm = Double.NaN; 
    CoordinatesPacker cp = new CoordinatesPacker(nRoots);
    
    for (int trial = 0; trial < 2; trial++) 
    {
      lognorm = Double.NEGATIVE_INFINITY;
      for (int i = 0; i < nRoots; i++)
        for (int j = 0; j < nRoots; j++) 
        {
          final int idx = cp.coord2int(i,j);
          if (i >= j) prs[idx] = Double.NEGATIVE_INFINITY;
          else
          {
            final boolean violation = SymmetricDiff.violatesOne(state.mergedClade(i,j), constraints);
            if (violation)
              prs[idx] = (trial == 0 ? Double.NEGATIVE_INFINITY : constraintViolationLogPenalty);
            if (trial == 0 && !violation)
            {
              double priorPostTerm = 0.0;
              if (usesPriorPost) // if not, it's uniform over allowed spans!
              {
                priorPostTerm = state.peekLogLikelihoodRatio(i,j,delta,0,0);
                lognorm = SloppyMath.logAdd(lognorm, priorPostTerm);
              }
              prs[idx] = priorPostTerm;
            }
            if (trial == 1)
            {
              double priorPostTerm = 0.0;
              if (usesPriorPost) // if not, it's uniform over allowed spans!
              {
                priorPostTerm = state.peekLogLikelihoodRatio(i,j,delta,0,0);
                lognorm = SloppyMath.logAdd(lognorm, priorPostTerm);
              }
              prs[idx] += priorPostTerm;
            }
          }
        }
      if (lognorm > Double.NEGATIVE_INFINITY) break;
    }
    
    NumUtils.expNormalize(prs);
    final int idx;
    try { idx = SampleUtils.sampleMultinomial(rand, prs); }
    catch (RuntimeException e) { return null; }
    final int [] pair = cp.int2coord(idx);
    PartialCoalescentState result = state.coalesce(pair[0],pair[1],delta,0,0);
    double w;
    if (usesPriorPost){ w = lognorm;                     }
    else              { w = result.logLikelihoodRatio(); } 
    return Pair.makePair(result, w);
  }
  public PartialCoalescentState getInitial() { return initial; }
  
  public int nIterationsLeft(PartialCoalescentState partialState)
  {
    return partialState.nIterationsLeft();
  }
}
