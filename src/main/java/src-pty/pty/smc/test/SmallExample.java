package pty.smc.test;
import java.io.*;
import java.util.*;

import pty.smc.ParticleFilter;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.ParticleKernel;

import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.exec.Execution;
import fig.prob.SampleUtils;
import gep.util.OutputManager;

import nuts.io.CSV;
import nuts.io.IO;
import nuts.tui.Table;
import nuts.util.CollUtils.*;
import nuts.util.CollUtils;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class SmallExample implements Runnable  
{
  @Option public File ranks;
  @Option public File qFwd;
  @Option public File qBwd;
  @Option public File gammaFile;
  
  @Option public double p = 0.5;
  
  @Option public int nReplicates = 1000;
  @Option public boolean noBack = false;
  
  private static ParticleFilter<Integer>pf = new ParticleFilter<Integer>();
  
  @Override
  public void run()
  {
    load();
    LogInfo.logs("Data loaded");
    // run n times
//    PrintWriter w = IOUtils.openOutHard(Execution.getFile("results.csv"));
//    w.println(CSV.header("replicate", "varDist"));
    OutputManager outMan = new OutputManager();
    for (int n = 0; n < nReplicates ; n++)
    {
      ParticleMapperProcessor<Integer, Integer> pp = ParticleMapperProcessor.saveParticlesProcessor();
      pf.sample(new SmallKernel(), pp);
      Counter<Integer> approx = pp.getCounter();
      approx.normalize();
      final double tvd = totalVariationDist(trueDistC, approx);
      final double iwv = impWeightVariance(approx);
      outMan.printWrite("variationalError", "replicate", n, "variationalError", tvd);
      outMan.printWrite("impWeightVariance", "replicate", n, "impWeightVariance", iwv);
//      LogInfo.logs("Approx:" + approx);
//      LogInfo.logs("Error:"+ tvd);
//      w.println(CSV.body(n, tvd));
//      if (n % 100 == 0)
//        w.flush();
    }
//    w.close();
  }
  
  private double impWeightVariance(Counter<Integer> approx)
  {
    final double mean = approx.totalCount() / ((double) approx.size());
    double sum = 0.0;
    for (double w : approx.entries.values())
      sum += (w - mean) * (w- mean);
    return sum / ((double) approx.size()); // - 1); avoid problem!
  }
  

  public static <T> double totalVariationDist(Counter<T> c1, Counter<T> c2)
  {
    double total = 0.0;
    for (T key : CollUtils.union(c1.keySet(), c2.keySet()))
      total += Math.abs(c1.getCount(key) - c2.getCount(key));
    return total * 0.5;
  }
  
  private int nRanks, nPStates;
  private double [][] fwd, bwd;
  private Map<Integer,Integer> i2r, r2i;
  private double [] gamma, trueDist;
  private Counter<Integer> trueDistC;
  
  public void load()
  {
    i2r = map();
    r2i = map();
    int r = 0;
    for (String line : IO.i(ranks)) if (!line.isEmpty())
    {
      String [] fields = line.split("[,]");
      for (String _i : fields)
      {
        nPStates++;
        int i = Integer.parseInt(_i);
        i2r.put(i, r);
        r2i.put(r, i);
      }
      r++;
      nRanks = r;
    }
    fwd = loadTrans(qFwd, true);
    bwd = loadTrans(qBwd, false);
    gamma = new double[nPStates];
    trueDist = new double[nPStates];
    trueDistC = new Counter<Integer>();
    for (String line : IO.i(gammaFile)) if (!line.isEmpty())
    {
      String [] split = line.split("\\t");
      int index = Integer.parseInt(split[0]);
      double value = Double.parseDouble(split[1]);
      gamma[index] = value;
      if (i2r.get(index) == nRanks - 1)
      {
        trueDistC.setCount(index, value);
        trueDist[index] = value;
      }
    }
    NumUtils.normalize(trueDist);
    trueDistC.normalize();
    LogInfo.logs("Gamma:" + Arrays.toString(gamma));
    LogInfo.logs("True dist:" + Arrays.toString(trueDist));
  }
  
  public class SmallKernel implements ParticleKernel<Integer>
  {

    @Override
    public Pair<Integer, Double> next(Random rand, Integer current)
    {
      // sample from q+
      int next = SampleUtils.sampleMultinomial(rand, fwd[current]);
      // get value of q+, q-, gammas
      double qPlus = fwd[current][next];
      double qMinus= (noBack ? 1.0 : bwd[next][current]);
      double newGamma = gamma[next];
      double oldGamma = gamma[current];
      double ratio = newGamma * qMinus / oldGamma / qPlus;
      return Pair.makePair(next, Math.log(ratio));
    }

    @Override
    public int nIterationsLeft(Integer partialState)
    {
      int rank = i2r.get(partialState);
      return nRanks - rank - 1;
    }

    @Override
    public Integer getInitial()
    {
      return 0;
    }
    
  }

  private double [][] loadTrans(File qFwd, boolean isFwd)
  {
    double[][] trans = new double[nPStates][nPStates];
    for (String line : IO.i(qFwd)) if (!line.isEmpty())
    {
      String [] fields = line.split("\\t");
      int f1 = Integer.parseInt(fields[0]);
      int f2 = Integer.parseInt(fields[1]);
      double value = parse(fields[2]);
      trans[isFwd ? f1 : f2][isFwd ? f2 : f1] = value;
    }
    NumUtils.normalizeEachRow(trans);
    LogInfo.logs((isFwd ? "Fwd" : "Bwd") + ":\n" + Table.fromMatrix(trans, true, true, false));
    return trans;
  }
  
  private double parse(String string)
  {
    if (string.equals("p")) return p;
    if (string.equals("q")) return 1.0 - p;
    return Double.parseDouble(string);
  }

  public static void main(String[] args)
  {
    IO.run(args, new SmallExample(), "pf", pf);
  }


}
