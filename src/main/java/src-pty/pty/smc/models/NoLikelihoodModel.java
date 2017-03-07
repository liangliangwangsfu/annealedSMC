package pty.smc.models;
import java.io.*;
import java.util.*;

import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.Taxon;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.Train;
import pty.eval.SymmetricDiff;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.ParticleFilter.ForkedProcessor;

/**
 * Use this to sample from the prior using the particle filter
 * @author bouchard
 */
public class NoLikelihoodModel implements LikelihoodModelCalculator
{

  public double logLikelihood() { return 0.0; }
  public double extendLogLikelihood(double delta)
  {
    return 0;
  }
  
  // sample from prior
  public static void main(String [] args)
  {
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    Execution.run(args, new SampleFromPrior(),
        "pcs", PartialCoalescentState.class,
        "pf", SampleFromPrior.pf);
  }
  
  public static class SampleFromPrior implements Runnable
  {
    @Option public ArrayList<String> leaves = new ArrayList<String>();
    public static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
    
    public void run() { sample(); }

    public RootedTree sample()
    {
      PartialCoalescentState init = PartialCoalescentState.initState(leaves);
      
      ParticleKernel<PartialCoalescentState> ppk = new PriorPriorKernel(init);
      
      ParticleFilter.ParticleMapperProcessor<PartialCoalescentState,PartialCoalescentState> mapd = 
           ParticleFilter.ParticleMapperProcessor.saveCoalescentParticlesProcessor();

      pf.sample ( ppk , mapd);
      PartialCoalescentState pcs = mapd.map(), pcs2 = mapd.sample(pf.rand);
      
      RootedTree map = pcs.getFullCoalescentState();
      RootedTree sampled =pcs2.getFullCoalescentState();
      Train.outputTree (map.topology(), "map-tree", map.branchLengths());
      LogInfo.logs("Log likelihood of map-tree:" + pcs.logLikelihood());
      Train.outputTree (sampled.topology(), "sampled-tree", sampled.branchLengths());
      LogInfo.logs("Log likelihood of map-tree:" + pcs2.logLikelihood());
      return pcs.getFullCoalescentState();
    }
    
  }

  public LikelihoodModelCalculator combine(LikelihoodModelCalculator node1,
      LikelihoodModelCalculator node2, double delta1, double delta2, boolean dnbc)
  {
    return this;
  }
  public boolean isReversible()
  {
    return true;
  }
  public double peekCoalescedLogLikelihood(LikelihoodModelCalculator node1,
      LikelihoodModelCalculator node2, double delta1, double delta2)
  {
    return 0;
  }
}
