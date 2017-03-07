package pty.smc.test;
import fig.basic.NumUtils;
import fig.basic.Pair;
import fig.prob.SampleUtils;
import hmm.Param;
import hmm.ParamUtils;

import java.io.*;
import java.util.*;

import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.test.TestParticleNormalization.HMMPState;

import nuts.math.Sampling;
import nuts.math.TabularGMFct;
import nuts.math.TreeSumProd;
import nuts.math.TreeSumProd.HmmAdaptor;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class TestConditionalSMC
{
  
  public static class HMMConditionalParticleKernel implements ParticleKernel<HMMPState>
  {
    public final Param params;
    public final int [] obs;
    
    public int left, right;
    public final int [] cStates;

    public HMMConditionalParticleKernel(Param params, int[] obs, int[] cStates)
    {
      this.params = params;
      this.obs = obs;
      this.cStates = cStates;
    }
    public HMMPState getInitial() { return new HMMPState(left-1, (right == 0 ? -1 : cStates[right-1]), null); }
    public int nIterationsLeft(HMMPState partialState)
    {
      return right - partialState.t - 1;
    }
    public Pair<HMMPState, Double> next(Random rand, HMMPState current)
    {
      if (current.t >= right) throw new RuntimeException();
      int nxt;
      if (current.t == right - 1 && right < cStates.length)
      {
        // need to take into account flanking symbol
        int rightSym = cStates[right];
        double [] unormPrs = new double[params.nStates()];
        for (int s = 0; s < params.nStates(); s++)
          unormPrs[s] = // 2 factors! : obs included in weights instead!
            params.transMtx.p(current.state, s) *
            params.transMtx.p(s, rightSym);
        NumUtils.normalize(unormPrs);
        nxt = SampleUtils.sampleMultinomial(rand, unormPrs);
      }
      else if (current.t >= 0)
        // sample according to the dynamics
        nxt = params.transMtx.nextState(current.state, rand);
      else
        nxt = params.initVec.nextState(rand);
      // weight is the observation:
      double w = Math.log(params.emiMtx.p(nxt,obs[current.t+1]));
      return Pair.makePair(new HMMPState(current.t+1,nxt, current),w);
    }
  }

  /**
   * @param args
   */
  public static void main(String[] args)
  {
    Random rand = new Random(1);
    // generate a random set of params
    Param p = ParamUtils.randomUniParam(rand,2,2);
//    System.out.println("Param:\n"+p);
    int [] obs = new int[20];
//    System.out.println("Observation:" + Arrays.toString(obs));
    
    // exact
    ArrayList<Integer> obsList = new ArrayList<Integer>();
    for (int o : obs) obsList.add(o);
    HmmAdaptor adapt = new HmmAdaptor(p,obsList);
    TreeSumProd<Integer> tsp = new TreeSumProd<Integer>(adapt);
//    System.out.println("Exact="+tsp.logZ());
    TabularGMFct<Integer> exactMoments = tsp.moments();
//    System.out.println("Exact marg:\n" + tsp.moments());
    
    // pf
//    smcApprox(10000, p, obs, rand, exactMoments, false);
//    System.out.println("----");
    for (int i = 10; i < 1000000000; i*=10)
    {
      System.out.println("I=" + i);
      smcApprox(i, p, obs, rand, exactMoments, true);
      smcApprox(i, p, obs, rand, exactMoments, false);
      System.out.println();
    }
    

  }

  private static void smcApprox(
      int maxIter, 
      Param p, 
      int[] obs, 
      Random rand, 
      TabularGMFct<Integer> exactMoments,
      boolean breakIt)
  {
    // init at random the unobserved
    int[] cStates = new int[obs.length];
    for (int i = 0; i < obs.length; i++)
      cStates[i] = (rand.nextInt(p.nObs()));
    
    ParticleFilter<HMMPState> pf = new ParticleFilter<HMMPState>();
    pf.N = 10;
    
//    final int maxIter = 100000;
    final int winSize = 2;//obs.length/2;
    
    Counter<String> stats = new Counter<String>();
    
    HMMConditionalParticleKernel kernel = new HMMConditionalParticleKernel(p, obs, cStates);
    for (int iter = 0; iter < maxIter; iter++)
    {
      kernel.left = rand.nextInt(obs.length - winSize + 1);
      kernel.right = kernel.left + winSize;
//      if (kernel.right >= obs.length)
//        throw new RuntimeException();
      // get current cStates and its weights
      List<HMMPState> currentCStates = list();
      double [] uw = new double[winSize];
      int _cur = 0;
      HMMPState prev = new HMMPState(kernel.left-1, -1, null);
      for (int t = kernel.left; t < kernel.right; t++)
      {
        int curCState = cStates[t];
        HMMPState curHMMPS = new HMMPState(t, curCState, prev);
        currentCStates.add(curHMMPS);
        uw[_cur] = p.emiMtx.p(curCState, obs[t]);
        _cur++;
        prev = curHMMPS;
      }
      // call PF
      if (!breakIt)
        pf.setConditional(currentCStates, uw);
      StoreProcessor<HMMPState> pro = new StoreProcessor<HMMPState>();
      pf.sample(kernel, pro);
      // sample one path (check if it's a new one just for stat purpose)
      int index = Sampling.sample(rand, pro.ws);
      HMMPState sampled = pro.particles.get(index);
      // set values in cState
      while (sampled.ancestor != null)
      {
        cStates[sampled.t] = sampled.state;
        sampled = sampled.ancestor;
      }
      // update suff stats
      for (int t = 0; t < obs.length; t++)
        stats.incrementCount("Z_" + t + "=" + cStates[t], 1.0);
    }
    
    double error = 0.0;
    for (int t = 0; t < obs.length; t++)
      for (int s = 0; s < p.nStates(); s++)
      {
//        System.out.println("Exact: " + exactMoments.get(t,s));
//        System.out.println("Approx:" + stats.getCount("Z_" + t + "=" + s)/maxIter);
        final double delta = Math.abs(exactMoments.get(t,s) - stats.getCount("Z_" + t + "=" + s)/maxIter);
        error += delta;
//        System.out.println();
      }
    System.out.println("Break=" + breakIt + ",error=" + error);
  }

}
