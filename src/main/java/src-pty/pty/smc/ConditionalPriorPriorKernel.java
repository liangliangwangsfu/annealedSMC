package pty.smc;
import java.io.*;
import java.util.*;

import pty.eval.SymmetricDiff;
import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.CoordinatesPacker;

import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.prob.SampleUtils;
import goblin.Taxon;

import pty.smc.MapLeaves;

/**
 */
public class ConditionalPriorPriorKernel implements ParticleKernel<PartialCoalescentState>
{
  @Option public static boolean usesPriorPost = false;
  private final PartialCoalescentState initial;
  private final Set<Set<Taxon>> conditionalClades;
  private double agreementWeight;
  private MapLeaves ml;
  public ConditionalPriorPriorKernel(
      PartialCoalescentState initial, 
      Set<Set<Taxon>> conditionalClades, 
      MapLeaves ml, 
      double agreementWeight) { 
	  this.agreementWeight = agreementWeight;
	  this.conditionalClades = ml.mapClades(conditionalClades);
	  this.ml = ml;
	  this.initial = initial;
	  if (!initial.isClock())
	    throw new RuntimeException("Not yet supported");
  }
  public Pair<PartialCoalescentState, Double> next (Random rand, PartialCoalescentState state) {
    int nRoots = state.nRoots();
    final double delta = Sampling.sampleExponential (rand, 1.0/PriorPriorKernel.nChoose2(nRoots));
    double [] prs = new double[nRoots*nRoots];
//    PartialCoalescentState [] nxtStates = (usesPriorPost ? new PartialCoalescentState[nRoots*nRoots] : null);
    double lognorm = (usesPriorPost ? Double.NEGATIVE_INFINITY : Double.NaN);
    CoordinatesPacker cp = new CoordinatesPacker(nRoots);
    Set<Set<Taxon>> allCladesInState;// = ml.filterClades(state.allClades());
    if (state.mapped !=  null) allCladesInState = state.mapped;
    else                       allCladesInState = ml.filterClades(state.allClades());
    for (int i = 0; i < nRoots; i++)
      for (int j = 0; j < nRoots; j++)
      {
        final int idx = cp.coord2int(i,j);
        if (i >= j) prs[idx] = Double.NEGATIVE_INFINITY;
        else
        {
          double priorPostTerm = 0.0;
          if (usesPriorPost)
          {
//            PartialCoalescentState pcs = state.coalesce(i,j, delta,0,0);
//            nxtStates[idx] = pcs;
            priorPostTerm = state.peekLogLikelihoodRatio(i,j,delta,0,0); //pcs.logLikelihoodRatio();
            lognorm = SloppyMath.logAdd(lognorm, priorPostTerm);
          }
          prs[idx] =  priorPostTerm - agreementWeight * SymmetricDiff.deltaSymmetricDiff(state, allCladesInState, i, j, conditionalClades, ml);
        }
      }
    NumUtils.expNormalize(prs);
    final int idx = SampleUtils.sampleMultinomial(rand, prs);
    final int [] pair = cp.int2coord(idx);
    PartialCoalescentState result; double w;
    result = state.coalesce(pair[0],pair[1],delta,0,0);
    if (usesPriorPost) w = lognorm;                     
    else               w = result.logLikelihoodRatio(); 
    
    // optimization
    Set<Set<Taxon>>  newAllClades = new HashSet<Set<Taxon>>(allCladesInState);
    newAllClades.add(ml.restrict(state.mergedClade(pair[0],pair[1])));
    result.mapped = newAllClades;
    
    return Pair.makePair(result, w);
  }
  public PartialCoalescentState getInitial() { return initial; }
  public int nIterationsLeft(PartialCoalescentState partialState)
  {
    return partialState.nIterationsLeft();
  }
}
