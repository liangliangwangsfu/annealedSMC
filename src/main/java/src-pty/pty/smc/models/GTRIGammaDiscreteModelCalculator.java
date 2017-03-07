/**
 * 
 */
package pty.smc.models;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import goblin.Taxon;
import ma.MSAPoset;
import ma.RateMatrixLoader;
import ma.SequenceType;
import ma.MSAPoset.Column;
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
 * 
 * @author L. Wang
 * 
 */
public class GTRIGammaDiscreteModelCalculator implements LikelihoodModelCalculator 
{
  @Option public static boolean allowDiscreteModelCalculator = false;  // set to true if you are sure you want to use the slow version of this
  public static Map<Taxon,GTRIGammaDiscreteModelCalculator> getInit(File alignmentFile, SequenceType sequenceType)
  {
    if (true)
      throw new RuntimeException("Use Dataset.java instead!");
    LogInfo.track("Loading data");
    // create the appropriate  encodings
    Indexer<Character> indexer = sequenceType.getEncodings().nonGapCharactersIndexer();

    // read alignment
    MSAPoset align = MSAPoset.parseAlnOrMsfFormats(alignmentFile);
    LogInfo.logs("Alignment: " + align.nTaxa() + " taxa and " + align.columns().size() + " sites");
    // for each full column, create a site
    List<Map<Taxon,Integer>> indices = CollUtils.list();
    for (Column c : align.columns())
      if (c.getPoints().size() == align.nTaxa())
      {
        Map<Taxon,Integer> current = CollUtils.map();
        Character lastChar = null;
        try
        {
          for (Taxon t : c.getPoints().keySet())
          {
            lastChar = align.charAt(c, t);
            current.put(t, indexer.o2i(lastChar));
          }
          indices.add(current);
        }
        catch (Exception e)
        {
          LogInfo.warning("Probably an unknown character:" + lastChar);
        }
      }
    
    // create CTMC
    int nSites = indices.size(); 
    LogInfo.logs("Kept " + nSites + " gap-free sites");
    CTMC ctmc = null;
    if (sequenceType == SequenceType.RNA || 
        sequenceType == SequenceType.DNA)
      ctmc = SimpleCTMC.dnaCTMC(nSites);
    else if (sequenceType == SequenceType.PROTEIN)
      ctmc = SimpleCTMC.proteinCTMC(nSites);
    else
      throw new RuntimeException();
    LogInfo.logs("" + sequenceType + " rate matrix loaded");
    
    // create the discreteModelCalculators from this
    Map<Taxon, GTRIGammaDiscreteModelCalculator> result = CollUtils.map();
    for (Taxon t : align.taxa())
    {
      double [][] init = new double[nSites][ctmc.nCharacter(0)];
      for (int s = 0; s < nSites; s++)
        init[s][indices.get(s).get(t)] = 1.0;
      result.put(t, observation(ctmc, init));
    }
    
    LogInfo.end_track();
    
    return result;
  }
  
  public LikelihoodModelCalculator combine(
      LikelihoodModelCalculator node1, LikelihoodModelCalculator node2, 
      double v1, double v2, boolean avoidBuildCache)
  {
    return (LikelihoodModelCalculator) calculate(node1, node2, v1, v2, false);
  }
  
  /// Can ignore stuff under there
  
  public static GTRIGammaDiscreteModelCalculator observation(CTMC ctmc, double [][] initCache, double alpha, int nGammaCat)
  {    
    initCache = Dataset.DatasetUtils.log(initCache);
    return new GTRIGammaDiscreteModelCalculator(initCache, ctmc,alpha, nGammaCat);
  }
  

  
  public double[] getCacheCopy(int site)
  {
    double [] result = new double[ctmc.nCharacter(site)];
    for (int i = 0; i < cache[site].length; i++)
      result[i] = cache[site][i];
    return result;
  }
  /*
   * site -> character (log scale)
   */
  private final ArrayList<double [][]> cache; 
//  private final boolean [] isMissing;
  private final double logLikelihood;
  public final CTMC ctmc;
  
	private final double alpha; 
	private final int nGammaCat; 
	private final double[] rates;

	
//  public boolean isMissing(int site) { return MathUtils.close(ctmc.nCharacter(site), SloppyMath.logAdd(cache[site])); }
  /**
   * @param cache
   * @param ctmc
   * @param logLikelihood
   */
  public GTRIGammaDiscreteModelCalculator(ArrayList<double[][]> cache, CTMC ctmc, double logLikelihood, double alpha, int nGammaCat)
  {
    if (!allowDiscreteModelCalculator)
      throw new RuntimeException("Should use FastDiscreteModelCalculator");
    this.cache = cache;
    this.ctmc = ctmc;
    this.logLikelihood = logLikelihood;
    this.alpha=alpha;
    this.nGammaCat=nGammaCat;
    this.rates=CTMC.GTRIGammaCTMC.calculateCategoryRates(nGammaCat, alpha, 0);
//    this.isMissing = isMissing;
  }
  /**
   * @param cache
   * @param ctmc
   */
  private GTRIGammaDiscreteModelCalculator(ArrayList<double[][]> cache, CTMC ctmc,double alpha, int nGammaCat) // call this if you are initializing !
  {
    if (!allowDiscreteModelCalculator)
      throw new RuntimeException("Should use FastDiscreteModelCalculator");
    this.cache = cache;
    this.ctmc = ctmc;
    this.alpha=alpha;
    this.nGammaCat=nGammaCat;
    this.rates=CTMC.GTRIGammaCTMC.calculateCategoryRates(nGammaCat, alpha, 0);
    this.logLikelihood =   extendLogLikelihood(0.0);
  }
  
  
  
  
  private static int getMAP(double[] cache)
  {
    int argmax = -1;
    double max = Double.NEGATIVE_INFINITY;
    for (int c =0; c < cache.length; c++)
    {
      double cval = cache[c];
      if (cval > max)
      {
        max = cval;
        argmax = c;
      }
    }
    return argmax;
  }

  

    
  public double extendLogLikelihood(double delta)
  {
    double logLikelihood = 0.0;
    double [][] pr = (ctmc.isSiteTied() ? ctmc.getTransitionPr(0,delta) : null);
    double [] initD= (ctmc.isSiteTied() ? ctmc.getInitialDistribution(0) : null);
    for (int s = 0; s < ctmc.nSites(); s++)
    {
      initD=(ctmc.isSiteTied() ? initD : ctmc.getInitialDistribution(s));
      pr = (ctmc.isSiteTied() ? pr : ctmc.getTransitionPr(s,delta));
      final double [] array = new double[ctmc.nCharacter(s)]; 
      for (int y = 0; y < ctmc.nCharacter(s); y++)
      {
        double sum = 0.0;
        for (int x = 0; x < ctmc.nCharacter(s); x++)
          //sum += pr[x][y] + initD[x];       //wrong
          sum += pr[x][y]*initD[x];         //correct
        array[y] = cache[s][y] + Math.log(sum);
      }
      final double siteLogLikelihood = SloppyMath.logAdd(array);
      if (isSiteLogLikelihoodValid(siteLogLikelihood)) 
        logLikelihood += siteLogLikelihood;
    }
    return logLikelihood;
  }


  
  /**
   * When combining two trees in a coalescent event,
   * return a fresh LikelihoodModelCalculator that
   * sits at the given height
   * @param node1
   * @param node2
   * @param height
   * @param isRoot Is the newly formed node the root of the tree?
   * @return
   */
  public Object calculate(
      final LikelihoodModelCalculator node1, final LikelihoodModelCalculator node2, 
      final double v1, final double v2,
      final boolean isPeek)
  {
    // resulting backward scores are indexed: site -> character, in log scale
    final double [][] result = (isPeek ? null :Dataset.DatasetUtils.createObsArray(ctmc));
//    boolean [] isMissing = new boolean[ctmc.nSites()];
    final GTRIGammaDiscreteModelCalculator 
      cache1 = (GTRIGammaDiscreteModelCalculator) (node1),
      cache2 = (GTRIGammaDiscreteModelCalculator) (node2);
    double logLikelihood = 0.0;
    // treat each indep. site
    final boolean isSiteTied = ctmc.isSiteTied();
    double [][] pr1 = (isSiteTied ? ctmc.getTransitionPr(0,v1) : null),
                pr2 = (isSiteTied ? ctmc.getTransitionPr(0,v2) : null);
    double [] initD = (isSiteTied ? ctmc.getInitialDistribution(0) : null);
//    LogInfo.logsForce("dmc prs1:" + Arrays.deepToString(pr1));
//    LogInfo.logsForce("dmc prs2:" + Arrays.deepToString(pr2));
    final int ns = ctmc.nSites();
    final double [] genericWorkingArray = (isSiteTied ?   new double[ctmc.nCharacter(0)] : null);
    for (int s = 0; s < ns; s++)
    {
      // compile transition probabilities
      // branch length is computed by (height of new node) - (height of child node)
      if (!isSiteTied)
      {
        pr1 = ctmc.getTransitionPr(s,v1);
        pr2 = ctmc.getTransitionPr(s,v2);
        initD=ctmc.getInitialDistribution(s);
      }
      // this temporary array is used for logAdd 
      final double [] array = (isSiteTied ? genericWorkingArray : new double[ctmc.nCharacter(s)]); 
      // v iterates over characters at the newly formed node
      double siteLogLikelihood = Double.NEGATIVE_INFINITY;
      final int ncs = ctmc.nCharacter(s);
      for (int v = 0; v < ncs; v++)
      {
        // take into account the prior over values of the root if we are at the root
        double prod = 0.0; //(isRoot ? Math.log(initD[v]) : 0.0);
        // take into account both left and right branches
        prod += sum(s, pr1[v], cache1, array); // log scale!
        prod += sum(s, pr2[v], cache2, array);
        if (!isPeek) result[s][v] = prod;
        siteLogLikelihood = SloppyMath.logAdd(siteLogLikelihood, prod + Math.log(initD[v]));
      }
      // not missing/observed with uncertainty if the siteLikelihood is in [0,1] 
      // (in logscale, allowing for slight numerical imprecision)
      if (isSiteLogLikelihoodValid(siteLogLikelihood)) 
        logLikelihood += siteLogLikelihood;
    }
    if (isPeek) return logLikelihood;
    else        return new GTRIGammaDiscreteModelCalculator(result, ctmc, logLikelihood);
  }


  private boolean isSiteLogLikelihoodValid(double number)
  {
    return number <= 0.000001;
  }
  /*
   * Computes \sum_w p(?->w) * cache(w)
   */
  private double sum(final int site, final double [] prs, final GTRIGammaDiscreteModelCalculator cacheAtLeaf, final double [] workingArray)
  {
    final double [] cal = cacheAtLeaf.cache[site];
    final int ncs = ctmc.nCharacter(site);
    for (int w = 0; w < ncs ; w++) 
      workingArray[w] = Math.log(prs[w]) + cal[w];
    return SloppyMath.logAdd(workingArray);
  }

  public double logLikelihood()
  {
    return logLikelihood;
  }

  public boolean isReversible()
  {
    return false;
  }

  public double peekCoalescedLogLikelihood(LikelihoodModelCalculator node1,
      LikelihoodModelCalculator node2, double delta1, double delta2)
  {
    return (Double) calculate(node1, node2, delta1, delta2, true);
  }



  
}