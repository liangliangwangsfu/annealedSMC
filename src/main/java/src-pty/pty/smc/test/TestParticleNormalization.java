package pty.smc.test;
import hmm.Param;
import hmm.ParamUtils;

import java.io.*;
import java.util.*;

import nuts.math.TreeSumProd;
import nuts.math.TreeSumProd.HmmAdaptor;

import fig.basic.Pair;
import fig.prob.SampleUtils;

import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.ParticleFilter.ResamplingStrategy;

public class TestParticleNormalization
{
  public static class HMMPState
  {
    public final int t;
    public final  int state;
    public final HMMPState ancestor;
    public HMMPState(int t, int state, HMMPState ancestor)
    {
      this.t = t;
      this.state = state;
      this.ancestor = ancestor;
    }
  }
  public static class HMMParticleKernel implements ParticleKernel<HMMPState>
  {
    public final Param params;
    public final int [] obs;
    public final int T;
    public HMMParticleKernel(Param params, int[] obs)
    {
      this.params = params;
      this.obs = obs;
      this.T = obs.length;
    }
    public HMMPState getInitial() { return new HMMPState(-1, 0, null); }
    public int nIterationsLeft(HMMPState partialState)
    {
      return T - partialState.t -1;
    }
    public Pair<HMMPState, Double> next(Random rand, HMMPState current)
    {
      if (current.t >= T) throw new RuntimeException();
      int nxt;
      if (current.t >= 0)
        // sample according to the dynamics
        nxt = params.transMtx.nextState(current.state, rand);
      else
        nxt = params.initVec.nextState(rand);
      // weight is the observation:
      double w = Math.log(params.emiMtx.p(nxt,obs[current.t+1]));
      return Pair.makePair(new HMMPState(current.t+1,nxt, current),w);
    }
  }
  public static void main(String [] args)
  {
    Random rand = new Random(1);
    // generate a random set of params
    Param p = ParamUtils.randomUniParam(rand,10,10);
    System.out.println("Param:\n"+p);
    int [] obs = new int[20];
    System.out.println("Observation:" + Arrays.toString(obs));
    // compute normalization using particule filter
    HMMParticleKernel pk = new HMMParticleKernel(p, obs);
    ParticleProcessor voidPro = new ParticleFilter.DoNothingProcessor();
    ParticleFilter pf = new ParticleFilter();
    pf.N = 100;
//    pf.resamplingStrategy = ResamplingStrategy.ALWAYS;
    pf.resamplingStrategy = ResamplingStrategy.ALWAYS;
    pf.nThreads=8;
    pf.resampleLastRound = false;
    pf.sample(pk, voidPro);
    System.out.println("Approx=" + pf.estimateNormalizer());
    // exact
    ArrayList obsList = new ArrayList();
    for (int o : obs) obsList.add(o);
    HmmAdaptor adapt = new HmmAdaptor(p,obsList);
    TreeSumProd tsp = new TreeSumProd(adapt);
    System.out.println("Exact="+tsp.logZ());
  }
}
