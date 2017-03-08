package pty.smc.test;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.Taxon;

import java.io.*;
import java.util.*;
import java.util.concurrent.CountDownLatch;

import nuts.util.CollUtils;
import nuts.util.Counter;

import pty.eval.SymmetricDiff;
import pty.io.HGDPDataset;
import pty.smc.ConstrainedKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.test.SymmetryTest.MeanHeight;
import pty.smc.test.TestBrownianModel.KernelType;

public class TestKernelsAgreement implements Runnable
{
  @Option public static double firstKernelNIterationFactor = 1.0;
  @Option public static KernelType kernel1 = KernelType.PRIOR_PRIOR;
  @Option public static KernelType kernel2 = KernelType.PRIOR_POST;
  @Option public static int baseN = 100;
  @Option public static int maxN = 100000;
  @Option public static int addToN = 0;
  @Option public static double multiplyToN = 2.0;
  @Option public static double variance = 1.5e-5;
  @Option public static boolean onlyOne = false;
  
  public static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
  public static void main(String [] args)
  {
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    Execution.run(args, new TestKernelsAgreement(),
        "pf", pf,
        "hgdp", HGDPDataset.class,
        "pcs", PartialCoalescentState.class,
        "ppk", PriorPriorKernel.class);
  }
  public void run()
  {
    PartialCoalescentState initGeneState = TestBrownianModel.initGeneState (variance);
    
    ParticleKernel<PartialCoalescentState> 
      ppk1 = kernel1.load(initGeneState, null),
      ppk2 = kernel2.load(initGeneState, null);
    
    for (int N = baseN; N < maxN; N = (int) (multiplyToN*N + addToN))
    {
      LogInfo.logs("N=" + N);
      ParticleFilter.ParticleMapperProcessor<PartialCoalescentState,Set<Set<Taxon>>> 
        mbr1 = SymmetricDiff.createCladeProcessor(),
        mbr2 = SymmetricDiff.createCladeProcessor();
      MeanHeight   
        mh1 = new MeanHeight(),
        mh2 = new MeanHeight();
      ForkedProcessor
        fp1 = new ForkedProcessor(mbr1, mh1),
        fp2 = new ForkedProcessor(mbr2, mh2);
      
      // sample first
      pf.N = (int) (N * firstKernelNIterationFactor);
      pf.sample(ppk1, fp1);
      if (!onlyOne)
      {
        // sample second
        pf.N = N;
        pf.sample(ppk2, fp2);
      }
      // check the difference in the suffstats
      LogInfo.logs(diff(mbr1.getCounter(), mbr2.getCounter()));
      LogInfo.logs("\nHeight\t" + mh1.ss.getSum() + "\t" + mh2.ss.getSum());
    }
  }
  
  @SuppressWarnings("unchecked")
  public static String diff(Counter c1, Counter c2)
  {
    Counter c = new Counter();
    for (Object o : CollUtils.union(c1.keySet(), c2.keySet()))
      c.setCount(o, Math.abs(c1.getCount(o) - c2.getCount(o)));
    StringBuilder result = new StringBuilder();
    for (Object o : c)
      result.append(o + "\t" + c1.getCount(o) + "\t" + c2.getCount(o) + "\t" + c.getCount(o) + "\n");
    return result.toString();
  }
}
