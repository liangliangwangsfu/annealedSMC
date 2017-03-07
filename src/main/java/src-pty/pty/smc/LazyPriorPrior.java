package pty.smc;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.Dataset.DatasetUtils;
import pty.learn.DiscreteBP;
import pty.mcmc.UnrootedTreeState;
import goblin.Taxon;
import java.io.*;
import java.util.*;

import ma.SequenceType;
import nuts.io.IO;
import nuts.math.GMFct;
import nuts.math.Sampling;
import nuts.math.TreeSumProd;
import nuts.util.Counter;
import nuts.util.MathUtils;

import ev.ex.InferenceExperiments;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;

import pty.smc.models.CTMC;

/**
 */
public class LazyPriorPrior implements ParticleKernel<LazyPCS>
{
  @Option public static boolean insureLastStepNoRootResampling = true;
  @Option public static LazyProposalType proposalType = LazyProposalType.COAL;
  @Option public static LazyProposalType priorType = LazyProposalType.COAL;
  private static final double uniformMax = 1.0;
//  @Option public static boolean useYule = false;
  
  public static enum LazyProposalType 
  { 
    UNIF {
      @Override public DeltaDistribution getInstance()
      {
        return new UniformDeltaDistribution();
      }
    },
    YULE {
      @Override public DeltaDistribution getInstance()
      {
        return new YuleDeltaDistribution();
      }
    },
    COAL {
      @Override public DeltaDistribution getInstance()
      {
        return new CoalescentDeltaDistribution();
      }
    };
//    , 
//    UNI_COAL_MIX {
//      @Override
//      public DeltaDistribution getProposal()
//      {
//        return null;
//      }
//    }, 
//    UNI_YULE_MIX {
//      @Override
//      public DeltaDistribution getProposal()
//      {
//        return null;
//      }
//    };
    public abstract DeltaDistribution getInstance();
  }
  
  private final DeltaDistribution q,prior;
  private final LazyPCS initial;
  public LazyPriorPrior(PartialCoalescentState initial) 
  { 
    this.q = proposalType.getInstance();
    this.prior = priorType.getInstance();
    this.initial = new LazyPCS(initial); 
  }
  public LazyPCS getInitial() { return initial; }
  public int nIterationsLeft(LazyPCS partialState)
  {
    return partialState.nIterationsLeft();
  }
  
  private static Dataset dataset = null;
  private static CTMC ctmc;
  
  public Pair<LazyPCS, Double> next(
      Random rand,
      LazyPCS _current)
  {
    PartialCoalescentState current = _current.getState();
    // 1- sample the exponential waiting time
    final double delta = q.sampleDelta(current, rand);//Sampling.sampleExponential(rand, 1.0/nChoose2(current.nRoots()));        
 //   System.out.println("rand: "+rand.toString()+"	rate: "+(initial.isRooted() ? 1.0 : 0.5)/nChoose2(current.nRoots()) +"	delta is: "+delta);
    
    // 2- sample a random pair (without replacement)
    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, current.nRoots(), 2);
    if (sampledIndices.size() != 2)
      throw new RuntimeException();
    
    final int i0 = sampledIndices.get(0),
              i1 = sampledIndices.get(1);
    
    double p;
//    if (current.nIterationsLeft() == 1 && insureLastStepNoRootResampling)
//    {
//      // horrible hack, should be fixed!
//      PartialCoalescentState merged = current.coalesce(i0, i1, delta, 0, 0);
//      RootedTree rt = merged.getFullCoalescentState();
//      UnrootedTree ut = UnrootedTree.fromRooted(rt);
//      if (dataset == null)
//      {
//        dataset = DatasetUtils.fromAlignment(InferenceExperiments.data, SequenceType.RNA);
//        ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
//      }
//      NonClockTreeState ncs = NonClockTreeState.initFastState(ut, dataset, ctmc);
//      p = ncs.logLikelihood();
//    }
//    else
    p = current.peekLogLikelihoodRatio(i0, i1, delta, 0,0);
    
    double logWeightUpdate = 
      + p 
      + prior.deltaLogDensity(current, delta) 
      - q.deltaLogDensity(current, delta);
    
    // 3- the weight update is simply equal to the ratio of the new likelihood score to the old one
    return Pair.makePair(new LazyPCS(current, delta, i0, i1), logWeightUpdate);
  }
  public static double nChoose2(double n) { return n*(n-1)/2; }

  public static interface DeltaDistribution
  {
    public double sampleDelta(PartialCoalescentState current, Random rand);
    public double deltaLogDensity(PartialCoalescentState current, double delta);
  }
  
//  public static class UniformMixProposal implements Proposal
//  {
//    public static double max = 0.5;
//    public static double mix = 0.5;
//    @Override
//    public double normalizedLogDensity(PartialCoalescentState current, double x)
//    {
//      IO.warnOnce("Using uniform mix activated!");
//      double result = mix * Math.exp(logPriorRatio(current, x));
//      if (x <= max)
//        result += (1.0 - mix)/max;
//      return Math.log(result);
//    }
//
//    @Override
//    public double propose(PartialCoalescentState current, Random rand)
//    {
//      if (rand.nextDouble() < mix)
//        return Sampling.sampleExponential(rand, 1.0/nChoose2(current.nRoots())); 
//      else
//        return Sampling.nextDouble(rand, 0.0, max);
//    }
//    
//  }
  
  
//  public static class UniformDeltaDistribution implements DeltaDistribution
//  {
//    @Override
//    public double deltaLogDensity(PartialCoalescentState current, double delta)
//    {
//      return 1.0/uniformMax;
//    }
//
//    @Override
//    public double sampleDelta(PartialCoalescentState current, Random rand)
//    {
//      return Sampling.nextDouble(rand, 0.0, uniformMax);
//    }
//  }
  
  public static class CoalescentDeltaDistribution implements DeltaDistribution
  {
    @Override
    public double deltaLogDensity(PartialCoalescentState current, double delta)
    {
      return Sampling.exponentialLogDensity(1.0/nChoose2(current.nRoots()), delta);
    }

    @Override
    public double sampleDelta(PartialCoalescentState current, Random rand)
    {
      return Sampling.sampleExponential(rand, 1.0/nChoose2(current.nRoots())); 
    }
  }
  
  public static class UniformDeltaDistribution implements DeltaDistribution
  {
    @Override
    public double deltaLogDensity(PartialCoalescentState current, double delta)
    {
      return Math.log(1.0/uniformMax);
    }

    @Override
    public double sampleDelta(PartialCoalescentState current, Random rand)
    {
      return Sampling.nextDouble(rand, 0.0, uniformMax);
    }
  }
  
  public static class YuleDeltaDistribution implements DeltaDistribution
  {
    @Override
    public double deltaLogDensity(PartialCoalescentState current, double delta)
    {
      return Sampling.exponentialLogDensity(1.0/current.nRoots(), delta);
    }

    @Override
    public double sampleDelta(PartialCoalescentState current, Random rand)
    {
      return Sampling.sampleExponential(rand, 1.0/current.nRoots()); 
    }
  }
  
  // ratio between current and next built with delta
//  public static double logPriorRatio(PartialCoalescentState current, double delta)
//  {
//    if (useYule)
//      return Sampling.exponentialLogDensity(1.0/current.nRoots(), delta);
//    else
//      return Sampling.exponentialLogDensity(1.0/nChoose2(current.nRoots()), delta);
//  }
}
