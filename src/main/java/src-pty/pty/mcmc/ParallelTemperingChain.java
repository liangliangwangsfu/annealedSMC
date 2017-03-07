package pty.mcmc;
import java.io.*;
import java.util.*;

import nuts.math.Sampling;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.EasyFormat;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.Parallelizer;
import fig.basic.UnorderedPair;
import fig.exec.Execution;
import goblin.Taxon;

public class ParallelTemperingChain
{
  private List<PhyloSampler> samplers;
  private List<SummaryStatistics> swapStats;  // size: samplers.size()-1
  private TemperingOptions options = _defaultTemperingOptions;
  private boolean initialized = false;
  private String outputPrefix = "";
  
  public int nChains() { return samplers.size(); }
  public PhyloSampler roomTempChain() { return samplers.get(0); }
  public PhyloSampler getChain(int i) { return samplers.get(i); }
  public TemperingOptions getOptions(){ return options; }
  
  public void init(PhyloSampler roomTemperatureSampler)
  {
    this.samplers = new ArrayList<PhyloSampler>(options.nChains);
    this.swapStats = new ArrayList<SummaryStatistics>(options.nChains);
    for (int i = 0; i < options.nChains; i++)
    {
      if (i == 0) this.samplers.add(roomTemperatureSampler);
      else        this.samplers.add(roomTemperatureSampler.createHeatedVersion(1.0+options.temperature * i * i));
      if (i < options.nChainMovesPerRoundPerChain - 1)
        this.swapStats.add(new SummaryStatistics());
      // create directory for trees to output
      File directory = getOutDir(i);
      directory.mkdir();
      samplers.get(i).setFileOutputPrefix(directory.getName() + "/");
    }
    initialized = true;
  }
  
  public File getOutDir(int i) { return new File(Execution.getFile(outputPrefix + "chain-" + i)); }
  
  @Override
  public String toString()
  {
    StringBuilder result = new StringBuilder("Output prefix: " + outputPrefix + "\n");
    for (int i = 0; i < nChains(); i++)
    {
      PhyloSampler current = samplers.get(i);
      result.append("#"  + i + ": " +
      		"T=" + EasyFormat.fmt2(current.getTemperature()) +", " + 
      		"LL=" + EasyFormat.fmt2(current.logLikelihood()) + ", " +
      		"maxLL=" + EasyFormat.fmt2(current.mleLogLikelihood()) + ", " +
      		"ratio=" + EasyFormat.fmt2(current.getMeanAcceptanceRatio()) + ", " +
      		"condAnn=" + EasyFormat.fmt2(current.getConditionalAnnealRatio()) + ", " +
      		"condFra=" + (current.isConditioning() ? EasyFormat.fmt2(current.getConditionalFraction()) : "n/a") + ", " +
      		"swapRatio=" + (i !=nChains() - 1 ? EasyFormat.fmt2(swapStats.get(i).getMean()) : "n/a") + ", " +
      		"details={" + current.detailedRatioToString()  + "}" 
      		+ "\n");
    }
    return result.toString();
  }

  public void sample()
  {
    long start = System.currentTimeMillis();
    if (!initialized) throw new RuntimeException();
    for (int iter = 0; iter < options.nRounds; iter++)
    {
      if ((iter+1) %options.printFreq == 0)
        log(iter);
      sampleChains();
      if (System.currentTimeMillis() - start > options.timeCutOff) 
      {
        closeFiles();
        return;
      }
      for (int swapIter = 0; swapIter < nSwapsPerRound(); swapIter++)
      {
        sampleSwap(options.rand);
        if (System.currentTimeMillis() - start > options.timeCutOff)
        {
          closeFiles();
          return;
        }
      }
    }
  }
  
  public void closeFiles()
  {
    for (PhyloSampler sampler : samplers)
      sampler.closeFile();
  }
  
  private void log(int iter)
  {
    LogInfo.track("Round " + iter + "/" + options.nRounds + 
        " (" + samplers.get(0).nIterations() + " " +
        		"chain iters, " + iter * nSwapsPerRound() + " swaps)");
//    Counter<UnorderedPair<Language,Language>> edgeCounts = 
//      roomTempChain().getEdgeUpdateCount();
//    LogInfo.logs("Coverage: " + edgeCounts.size() + "/"
//        + roomTempChain().getInitialState().getNonClockTree(). + " " + edgeCounts);
//    System.out.println(roomTempChain().getInitialState().getNonClockTree().nonTerminalEdges());
//    System.out.println("Symm diff:" + CollUtils.symmetricDifference(
//        new HashSet(roomTempChain().getInitialState().getNonClockTree().nonTerminalEdges()),
//        edgeCounts.keySet()));
    LogInfo.logs(this.toString());
    LogInfo.end_track();
  }

  public int nSwapsPerRound() { return options.nSwapsPerRoundPerChain * nChains(); }
  
  private void sampleSwap(Random rand)
  {
    if (nChains() == 1) return;
    final int i = rand.nextInt(nChains()-1);
    PhyloSampler.sampleSwap(samplers.get(i),samplers.get(i+1),rand,swapStats.get(i));
  }
  
  private void sampleChains()
  {
    Parallelizer<Integer> parallelizer = new Parallelizer<Integer>(options.nThreads);
    parallelizer.setPrimaryThread();
    parallelizer.process(CollUtils.ints(nChains()), new Parallelizer.Processor<Integer>() {
      public void process(Integer x, int _i, int _n, boolean log) {
        for (int i = 0; i < options.nChainMovesPerRoundPerChain; i++)
          samplers.get(x).sample();
      }
    });
  }
  
  public static TemperingOptions _defaultTemperingOptions = new TemperingOptions();
  public static class TemperingOptions
  {
    @Option public long timeCutOff = Long.MAX_VALUE;
    @Option public int nThreads = 1;
    @Option public int nChains = 4;
    @Option public double temperature = 0.2;
    @Option public int nRounds = 100000;
    @Option public int nSwapsPerRoundPerChain = 2;
    @Option public int nChainMovesPerRoundPerChain = 10;
    @Option public int printFreq = 100;
    @Option public Random rand = new Random(1);
  }
  /**
   * returns the directory that will contain unheated samples
   * @param outputPrefix
   * @return
   */
  public File setOutputPrefix(String outputPrefix)
  {
    this.outputPrefix = outputPrefix;
    return getOutDir(0);
  }
}
