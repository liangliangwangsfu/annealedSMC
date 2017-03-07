/**
 * 
 */
package pty.smc.models;
import java.io.File;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.lucene.util.OpenBitSet;

import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Pair;
import fig.prob.SampleUtils;
import goblin.Taxon;
import ma.MSAPoset;
import ma.RateMatrixLoader;
import ma.SequenceType;
import ma.MSAPoset.Column;
import nuts.io.IO;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils;
import nuts.util.CounterMap;
import nuts.util.Indexer;
import nuts.util.MathUtils;
import pepper.Encodings;
import pty.io.Dataset;
import pty.smc.PartialCoalescentState;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.models.CTMC.SimpleCTMC;

/**
 * Likelihood models based on finite state CTMC
 * Has methods to initialize a likelihood model based on observations,
 * to update the likelihoods for a coalescent event.
 * 
 * Uses scalings instead of logs
 * @author bouchard
 */
public final class FastDiscreteModelCalculator implements LikelihoodModelCalculator, Serializable 
{
  private static final long serialVersionUID = 1L;

  public LikelihoodModelCalculator combine(
      LikelihoodModelCalculator node1, LikelihoodModelCalculator node2, 
      double v1, double v2, boolean avoidBuildCache)
  {
    return (LikelihoodModelCalculator) calculate(node1, node2, v1, v2, false, avoidBuildCache);
  }
  
  public static FastDiscreteModelCalculator observation(CTMC ctmc, double [][] initCache)
  {    
    return new FastDiscreteModelCalculator(initCache, ctmc, false);
  }
  
//  public static FastDiscreteModelCalculator observationWithItemizedLogLikelihood(CTMC ctmc, double [][] initCache)
//  {    
//    double [] array = new double[ctmc.nSites()];
//    return new FastDiscreteModelCalculator(initCache, ctmc, false, array);
//  }
  
  
  public static FastDiscreteModelCalculator observation(CTMC ctmc, double [][] initCache, boolean resampleRoot)
  {    
    return new FastDiscreteModelCalculator(initCache, ctmc, resampleRoot);
  }
 
  /*
   * site -> character (std scale)
   */
  private final double [][] cache; 
  private final double logLikelihood;
//  public final double [] itemizedLogLikelihoods;
  public final CTMC ctmc;
//  private final boolean resampleRoot;
  private final OpenBitSet [] bitVectorVersion;
  
  
  
//  private final DiscreteModelCalculator tester;
  
  public double[][] getCache()
  {
    return cache;
  }

  // called when combining
  private FastDiscreteModelCalculator(double[][] cache, CTMC ctmc,
      double logLikelihood, boolean resampleRoot) //, double [] itemizedLogLikelihoods)//, DiscreteModelCalculator tester)
  {
    this.bitVectorVersion = null;
//    this.resampledVersion = null; //resampleToBitVector(cache);
//    this.tester = tester;
    this.cache = cache;
    this.ctmc = ctmc;
    this.logLikelihood = logLikelihood;
//    this.itemizedLogLikelihoods = itemizedLogLikelihoods;
//    this.resampleRoot = resampleRoot;
  }

  // called when initialized
  private FastDiscreteModelCalculator(double[][] cache, CTMC ctmc, boolean resampleRoot) //, double [] itemizedLogLikelihoods) // call this if you are initializing !
  {
    if (resampleRoot)
      throw new RuntimeException("Discontinued feature");
    
//    IO.warnOnce("Using instrumented FDM... change that!");
//    this.tester = DiscreteModelCalculator.observation(ctmc,cache);
    this.ctmc = ctmc;
//    this.itemizedLogLikelihoods = itemizedLogLikelihoods;
    this.logLikelihood = initLogLikelihood(cache, ctmc); //, itemizedLogLikelihoods);
    this.cache = initCache(cache, ctmc);
    this.bitVectorVersion = observationToBitVector(cache); //resampleToBitVector(this.cache);
//    this.resampleRoot = resampleRoot;
    
  }

  

  private static double[][] initCache(double[][] cache, CTMC ctmc)
  {
//    return cache;
    if (!ctmc.isSiteTied())
      throw new RuntimeException();
    double [] initD = ctmc.getInitialDistribution(0);
    double [][] result = new double[cache.length][];
    for (int i = 0; i < result.length; i++)
    {
      // this handles unobserved or partially observed r.v.s
      result[i] = new double[cache[i].length];
      for (int j = 0; j < result[i].length; j++)
        result[i][j] = initD[j] * cache[i][j];
      // Note: the line below is important, caused a bug without it
      NumUtils.normalize(result[i]);
    }
    return result;
  }


//  private static double initLogLikelihood(double[][] cache, CTMC ctmc) //, double [] itemized)
//  {
//    if (!ctmc.isSiteTied())
//      throw new RuntimeException();
//    
////    final boolean needItemized = itemized != null;
//    
//    double [] initD = ctmc.getInitialDistribution(0);
//    double result = 0.0;
//    for (int s = 0; s < cache.length; s++)
//      for (int c = 0; c < initD.length; c++)
//      {
//        final double cur = cache[s][c];
//        if (cur != 0.0 && cur != 1.0)      
//          throw new RuntimeException();
//        final double term = cur * Math.log(initD[c]);
////        if (needItemized)
////          itemized[s] += term;
//        result += term;
//      }
//    return result;
//  }
// 
  
  private static double initLogLikelihood(double[][] cache, CTMC ctmc) //, double [] itemized)
  {
    if (!ctmc.isSiteTied())
      throw new RuntimeException();
    
//    final boolean needItemized = itemized != null;
    
    double [] initD = ctmc.getInitialDistribution(0);
    double result = 0.0;
    for (int s = 0; s < cache.length; s++)
    {
    	double term=0; 
      for (int c = 0; c < initD.length; c++)
      {
        final double cur = cache[s][c];
        if (cur != 0.0 && cur != 1.0)      
          throw new RuntimeException();
             term+= cur * initD[c];
//        if (needItemized)
//          itemized[s] += term;        
      }
      result += Math.log(term);
    }
    return result;
  }
 
  
  
  public double extendLogLikelihood(double delta)
  {
    throw new RuntimeException();

  }
  
  
  public static double quickPeek(      
      final LikelihoodModelCalculator node1, final LikelihoodModelCalculator node2, 
      final double v1, final double v2)
  {
    final FastDiscreteModelCalculator 
    cache1 = (FastDiscreteModelCalculator) (node1),
    cache2 = (FastDiscreteModelCalculator) (node2);
    final double [][] pr = cache1.ctmc.getTransitionPr(0,v1+v2);
    final double [] initD = cache2.ctmc.getInitialDistribution(0);
    final double [] logInit = new double[initD.length];
    for (int i = 0; i < logInit.length; i++)
      logInit [i]= Math.log(initD[i]);

    final int ncs = cache1.ctmc.nCharacter(0);
    
    double sum = 0.0;
    for (int x = 0; x < ncs; x++)
    {
      final double cur = logInit[x];
      final double [] curAr = pr[x];
      for (int xPrime = 0; xPrime < ncs; xPrime++)
      {
        final long n = OpenBitSet.intersectionCount(cache1.bitVectorVersion[x], cache2.bitVectorVersion[xPrime]);
        sum += n * (cur + Math.log(curAr[xPrime]));
      }
    }
    return sum;
  }
//  
//  public static void test(final LikelihoodModelCalculator node1, final LikelihoodModelCalculator node2)
//  {
//    final int nPart = 10000;
//    double delta = 0.0001;
//    
//    final FastDiscreteModelCalculator 
//    cache1 = (FastDiscreteModelCalculator) (node1);
//    
//    if (cache1.resampledVersion == null)
//      return;
//    
//    long timeForExacts = 0, timeForApprox = 0;
//    
//    double [] exact =  new double[nPart], approx = new double[nPart];
//    
//    for (int i = 0; i < nPart; i++)
//    {
//      double cur = (i+1)*delta;
//      long curTime = System.currentTimeMillis();
//      double exactValue = ((Double) cache1.calculate(node1, node2, cur/2.0, cur/2.0, true, true)).doubleValue();
//      timeForExacts += (System.currentTimeMillis() - curTime);
//      curTime = System.currentTimeMillis();
//      double approxValue = quickApproximateCalculate(node1, node2, cur/2.0, cur/2.0);
//      timeForApprox += (System.currentTimeMillis() - curTime);
//      exact[i] = exactValue;
//      approx[i] = approxValue;
//    }
//    NumUtils.normalize(exact);
//    NumUtils.normalize(approx);
//    
//    System.out.println("Time for exact:" + timeForExacts);
//    System.out.println("Time for approx:" + timeForApprox);
//    
//    System.out.println("Exact\tApprox");
//    for (int i = 0 ; i < nPart; i++)
//      System.out.println("" + ((i+1)*delta) + "\t"+ exact[i] + "\t"+ approx[i]);
//  }
  
  /**
   * When combining two trees in a coalescent event,
   * return a fresh LikelihoodModelCalculator that
   * sits at the given height
   * @param node1
   * @param node2
   * @param height
   * @return
   */
  private final Object calculate(
      final LikelihoodModelCalculator node1, final LikelihoodModelCalculator node2, 
      final double v1, final double v2,
      final boolean isPeek, final boolean avoidBuildCache)
  {
    
    if (isPeek && !avoidBuildCache)
      throw new RuntimeException();
    
    // resulting backward scores are indexed: site -> character, in std scale
    final double [][] result = (avoidBuildCache ? null :Dataset.DatasetUtils.createObsArray(ctmc));
    final FastDiscreteModelCalculator 
      cache1 = (FastDiscreteModelCalculator) (node1),
      cache2 = (FastDiscreteModelCalculator) (node2);
    
//    final boolean needItemized = !avoidBuildCache && (cache1.itemizedLogLikelihoods != null);
//    if (needItemized != (cache2.itemizedLogLikelihoods != null))
//      throw new RuntimeException();
    
    if (isPeek && cache1.bitVectorVersion != null && cache2.bitVectorVersion != null)
      return quickPeek(node1, node2, v1, v2);
    
//    DiscreteModelCalculator dmc = (DiscreteModelCalculator) tester.combine(cache1.tester, cache2.tester, v1, v2, avoidBuildCache);
//    System.out.println("DMC: " + dmc.logLikelihood());
    
//    final double [] newItemized = needItemized ? new double[ctmc.nSites()] : null;
    
    final double [][] pr1 = ctmc.getTransitionPr(0,v1),
                pr2 = ctmc.getTransitionPr(0,v2);
    final double [] initD = ctmc.getInitialDistribution(0);
    final double [] invD = new double[initD.length];
    
    for (int i = 0; i < invD.length; i++)
      invD[i] = 1.0/initD[i];
    
//    {
//      System.out.println(Arrays.deepToString(pr1));
//      System.out.println(Arrays.deepToString(pr2));
//      System.out.println(Arrays.toString(initD));
//      System.out.println(Arrays.toString(invD));
//    }
    
    
    final int ns = ctmc.nSites();
    
    
    final int ncs = ctmc.nCharacter(0);
    double logNorm = 0.0, tempNorm = 1.0;
    for (int s = 0; s < ns; s++)
    {
      // compute the norm of result[s][x]
      double currentNorm = 0.0;
      // optionally, fill an unormalized version of cache[s][x]...
      final double [] c1 = cache1.cache[s], c2 = cache2.cache[s];

      for (int x = 0; x < ncs; x++)
      {
        final double [] cPr1 = pr1[x], cPr2 = pr2[x];
        double s1 = 0.0, s2 = 0.0;
        for (int xPrime = 0; xPrime < ncs; xPrime++)
        {
          final double cInvD = invD[xPrime];
          final double cc1 = c1[xPrime]; if (cc1 != 0.0) s1 += cPr1[xPrime] * cInvD * cc1;
          final double cc2 = c2[xPrime]; if (cc2 != 0.0) s2 += cPr2[xPrime] * cInvD * cc2;
        }
        final double currentProd = initD[x] * s1 * s2;
        if (!avoidBuildCache) result[s][x] = currentProd;
        currentNorm += currentProd;
      }
//      if (needItemized)
//        newItemized[s] = Math.log(currentNorm) + cache1.itemizedLogLikelihoods[s] + cache2.itemizedLogLikelihoods[s];
      if (!avoidBuildCache) 
        MathUtils.normalizeAndGetNorm(result[s]);
      tempNorm *= currentNorm;
      double abs = Math.abs(tempNorm);
      if (abs < 1e-100 || abs > 1e+100 || s == ns -1)
      {
//        {
//          System.out.println("\t"+ tempNorm);
//        }
        
        logNorm += Math.log(tempNorm);
        tempNorm = 1.0;
      }
    }
    // NB: exp(lognorm) is equal to P(y_{left+right})/P(y_left)/P(y_right)
    final double newLogLikelihood = logNorm + cache1.logLikelihood + cache2.logLikelihood;
    // hence, to get P(y_{left+right}), can take the product (sum the logs)
    // i.e.: incremental telescoping product
    
//    System.out.println("Ours: " + newLogLikelihood);
//    System.out.println();
//    MathUtils.checkClose(newLogLikelihood, dmc.logLikelihood());
    
//    System.out.println("cache1.logLikelihood" +cache1.logLikelihood + "cache2.logLikelihood"+ cache2.logLikelihood+"newLogLikelihood "+newLogLikelihood);
    if (isPeek) return newLogLikelihood;
    else        return //resampleRoot && !avoidBuildCache ? null : //resampleRoot(result, ctmc) :
                                      new FastDiscreteModelCalculator(result, ctmc, newLogLikelihood, false); 
//        dmc);
  }

//  private static FastDiscreteModelCalculator resampleRoot(double[][] message, CTMC ctmc)
//  {
//    Random rand = myRandom.get();
//    for (int site = 0; site < message.length; site++)
//    {
//      int idx = SampleUtils.sampleMultinomial(rand, message[site]);
//      for (int i = 0; i < message[site].length; i++)
//        if (i == idx)
//          message[site][i] = 1.0;
//        else
//          message[site][i] = 0.0;
//    }
//    return new FastDiscreteModelCalculator(message, ctmc, true);
//  }
//  
//  private static ThreadLocal<Random> myRandom = new ThreadLocal<Random>()
//  {
//    public Random initialValue() { return new Random(1); }
//  };
  
//  private static OpenBitSet [] resampleToBitVector(double [][] message)
//  {
//    if (message == null)
//      return null;
//    Random rand = myRandom.get();
//    final int nChars = message[0].length;
//    final int nSites = message.length;
//    OpenBitSet [] result = new OpenBitSet[nChars];
//    for (int i = 0; i < result.length; i++)
//      result[i] = new OpenBitSet(nSites);
//    for (int site = 0; site < message.length; site++)
//    {
//      int idx = SampleUtils.sampleMultinomial(rand, message[site]);
//      result[idx].fastSet(site);
//    }
//    return result;
//  }

  private static OpenBitSet [] observationToBitVector(double [][] message)
  {
    final int nChars = message[0].length;
    final int nSites = message.length;
    OpenBitSet [] result = new OpenBitSet[nChars];
    for (int i = 0; i < result.length; i++)
      result[i] = new OpenBitSet(nSites);
    for (int site = 0; site < message.length; site++)
    {
      int idx = findOne(message[site]);
      if (idx == -1)
        return null;
      result[idx].fastSet(site);
    }
    return result;
  }
  
  
  private static int findOne(double[] ds)
  {
    int result = -1;
    for (int i = 0; i < ds.length; i++)
    {
      final double d = ds[i];
      if (d == 1)
      {
        if (result != -1)
          return -1; // more than one == 1
        result = i;
      }
    }
    return result;
  }

  public double logLikelihood()
  {
    return logLikelihood;
  }

  public boolean isReversible()
  {
    throw new RuntimeException();
  }

  public double peekCoalescedLogLikelihood(LikelihoodModelCalculator node1,
      LikelihoodModelCalculator node2, double delta1, double delta2)
  {
    return (Double) calculate(node1, node2, delta1, delta2, true, true);
  }
}