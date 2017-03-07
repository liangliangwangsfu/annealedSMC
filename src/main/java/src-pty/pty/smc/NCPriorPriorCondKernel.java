//package pty.smc;
//import java.io.*;
//import java.util.*;
//
//import fig.basic.Option;
//import fig.basic.Pair;
//import nuts.math.Sampling;
//import nuts.util.CollUtils.*;
//import nuts.util.Counter;
//import static nuts.util.CollUtils.*;
//import static nuts.io.IO.*;
//import static nuts.util.MathUtils.*;
//
//
//public class NCPriorPriorCondKernel implements ParticleKernel<PartialCoalescentState>
//{
//  
//  private final PartialCoalescentState initial;
//  
//  public NCPriorPriorCondKernel(PartialCoalescentState initial) { this.initial = initial; }
//  
//  @Override
//  public PartialCoalescentState getInitial() { return initial; }
//  
//  @Override
//  public int nIterationsLeft(PartialCoalescentState partialState)
//  {
//    return partialState.nIterationsLeft();
//  }
//  
//  @Override
//  public Pair<PartialCoalescentState, Double> next(Random rand,
//      PartialCoalescentState current)
//  {
//    // 1-sample branch lengths
//    final double delta0 = Sampling.sampleExponential(rand, 1.0/NCPriorPriorKernel.deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0),
//                 delta1 = Sampling.sampleExponential(rand, 1.0/NCPriorPriorKernel.deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0);
//       
//     // 2- sample a random pair (without replacement)
//    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, current.nRoots(), 2);
//    final int i0 = sampledIndices.get(0),
//              i1 = sampledIndices.get(1);
//    
//    final PartialCoalescentState next = current.coalesce(i0, i1, 0.0, delta0, delta1);
//    
//    // 3 - compute LOG weight..
//    final double logW = next.nonClockLogWeight();
//       
//    return Pair.makePair(next, logW);
//  }
//}

/////


package pty.smc;
import java.io.*;
import java.util.*;

import pty.eval.SymmetricDiff;
import pty.smc.PartialCoalescentState.CoalescentNode;

import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.prob.SampleUtils;
import goblin.Taxon;
import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils.*;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.CoordinatesPacker;
import nuts.util.Counter;
import nuts.util.MathUtils;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class NCPriorPriorCondKernel implements ParticleKernel<PartialCoalescentState>
{
  
  private final PartialCoalescentState initial;
  private Set<Set<Taxon>> conditioningClades;
  public int getNMasks() { return conditioningClades.size(); }
  public void setConditioning(Set<Set<Taxon>> conditioningClades)
  {
    this.conditioningClades = conditioningClades;
  }
  
  public NCPriorPriorCondKernel(PartialCoalescentState initial) 
  { 
    this.conditioningClades = set();
    this.initial = initial; 
  }
  
  @Override
  public PartialCoalescentState getInitial() { return initial; }
  
  @Override
  public int nIterationsLeft(PartialCoalescentState partialState)
  {
    return partialState.nIterationsLeft();
  }
  
//  private static final ThreadLocal < int[] > clades
  
  @Override
  public Pair<PartialCoalescentState, Double> next(Random rand,
      PartialCoalescentState current)
  {
    // 1-sample branch lengths
    final double delta0 = Sampling.sampleExponential(rand, 1.0/NCPriorPriorKernel.deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0),
                 delta1 = Sampling.sampleExponential(rand, 1.0/NCPriorPriorKernel.deltaProposalRate)/(current.nRoots()==2 ? 2.0 : 1.0);
    
    int nPossibleMergeInT = 0;
    int nRoots = current.nRoots();
    double [] prs = new double[nRoots*nRoots];
    CoordinatesPacker cp = new CoordinatesPacker(nRoots);
    for (int i = 0; i < nRoots; i++)
      for (int j = i+1; j < nRoots; j++)
        if (!SymmetricDiff.violatesOne(CollUtils.union(current.rootedClade(i), current.rootedClade(j)), conditioningClades))
        {
          prs[cp.coord2int(i,j)] = 1.0;
          nPossibleMergeInT++;
        }
    NumUtils.normalize(prs);
    int idx = -1;
    try {
     idx = SampleUtils.sampleMultinomial(rand, prs);
    } catch (Exception e)
    {
      return Pair.makePair(null, Double.NEGATIVE_INFINITY);
    }
    final int [] pair = cp.int2coord(idx);
    final int i0 = pair[0], i1 = pair[1];
    
    final PartialCoalescentState next = current.coalesce(i0, i1, 0.0, delta0, delta1);
    
//    final int nNonTrivialRootsInTPlus1 = next.nNonTrivialRoots();
    // 3 - compute LOG weight..
    final double logW = nonClockLogWeight(next);
       
    return Pair.makePair(next, logW);
  }
  
  public double nonClockLogWeight(PartialCoalescentState next)
  {
    if (next.isClock()) throw new RuntimeException();
    double logSum = Double.NEGATIVE_INFINITY;
    for (Arbre<CoalescentNode>  root : next.roots)
      if (!root.isLeaf())
      {
        if (root.getChildren().size() != 2)
          throw new RuntimeException();
        final double logCur = + root.getChildren().get(0).getContents().likelihoodModelCache.logLikelihood() 
                              + root.getChildren().get(0).getContents().likelihoodModelCache.logLikelihood()  
                              - root.getContents().likelihoodModelCache.logLikelihood() 
                              - topoNorm(next, root);
        logSum = SloppyMath.logAdd(logSum, logCur);
      }
    return -logSum;
  }

  private double topoNorm(PartialCoalescentState next, Arbre<CoalescentNode> root)
  {
    // a prior, without constraints:
    int result = MathUtils.safeIntValue(MathUtils.nChoose2(next.nRoots() + 1));
    // remove values due to constraints
    if (!root.isLeaf())
    {
      result -= nViolations(next, root.getChildren().get(0), root);
      result -= nViolations(next, root.getChildren().get(1), root);
    }
    return result;
  }

  private int nViolations(PartialCoalescentState next,
      Arbre<CoalescentNode> child, Arbre<CoalescentNode> parent)
  {
    int nViol = 0;
    for (Arbre<CoalescentNode> cRoot : next.roots)
      if (cRoot != parent)
        if (SymmetricDiff.violates(child.getContents().rootedClade, cRoot.getContents().rootedClade))
          nViol++;
    return nViol;
  }
  
}