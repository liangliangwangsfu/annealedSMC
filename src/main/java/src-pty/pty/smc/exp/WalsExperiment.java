package pty.smc.exp;
import java.io.*;
import java.util.*;

import nuts.tui.Table;
import nuts.util.Arbre;
import nuts.util.Counter;
import nuts.util.MathUtils;

import pty.eval.Purity;
import pty.eval.SymmetricDiff;
import pty.io.WalsAnn;
import pty.io.WalsData;
import pty.io.WalsDataset;
import pty.learn.CTMCExpectations;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;

import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.exec.Execution;
import fig.prob.SampleUtils;
import goblin.Taxon;

public class WalsExperiment //implements Runnable
{
//  @Option public double initRate = 1;
//  @Option public int nEMIters = 20;
//  @Option public boolean generate = true;
//  @Option public int generateDepth = 2;
//  @Option public int generateLength = 50;
//  @Option public Random generateRand = new Random(1);
//  
//  private ParticleFilter<PartialCoalescentState> pf;
//  
//  public WalsExperiment(ParticleFilter<PartialCoalescentState> pf) { this.pf = pf; }
//  
//  public void run()
//  {
////    WalsDataset data = (generate ? null : WalsDataset.getPreprocessedCorpus());
////    // create observations
////    Map<Language, int[]> observations = 
////      ( generate ? generate(generatingParams(generateLength), generateDepth, generateRand) :
////                   data.toObservationArrays(-1)); 
////    // estimate stationary distribution
////    double [][] statDistn = statDistn(observations);
////    //  rate matrix
////    CTMC ctmc = (generate ? generatingParams(generateLength) : ctmc(statDistn, initRate));
////    for (int emIter = 0; emIter < nEMIters; emIter++)
////    {
////      LogInfo.track("EM Iteration " + (emIter+1) + "/" + nEMIters,true);
////      // init state
////      PartialCoalescentState initState = initState(observations, ctmc);
////      // prepare sampling machinery
////      MinBayesRiskDecoder<PartialCoalescentState, Set<Set<Language>>> mbr =
////        SymmetricDiff.getApproximateMBRProcessor();
////      CTMCExpectations.SuffStatProcessor ssp = new CTMCExpectations.SuffStatProcessor(ctmc);
////      ForkedProcessor<PartialCoalescentState> processors = new ForkedProcessor<PartialCoalescentState>();
////      processors.processors.add(mbr);
////      processors.processors.add(ssp);
////      PriorPriorKernel kernel = new PriorPriorKernel(initState);
////      // sample!
////      pf.bootstrapFilter(kernel, processors);
////      Arbre<Language> reconstruction = SymmetricDiff.clades2arbre(mbr.centroid());
////      LogInfo.logs("Reconstructed tree:\n" + reconstruction.deepToString());
////      // evaluate purity
////      if (!generate)
////      {
////        Map<Language,String> lang2genus = WalsDataset.langDB.genusMap(); //lang2genus();//data);
////        LogInfo.logs("Purity:" + Purity.purity(Arbre.arbre2Tree(reconstruction), lang2genus));
////      }
////      ctmc = (generate ? ssp.getTiedMLE() : ssp.getMLE());
////      if (generate)
////        LogInfo.logs("Param:\n" + ctmc.toString());
////      LogInfo.end_track();
////    }
//  }
//  
//  public static Map<Language,int[]>  generate(CTMC model, int depth, Random rand)
//  {
//    if (depth <= 1) throw new RuntimeException();
//    LogInfo.logs("Generating from:\n" + model.toString());
//    Map<Language,int[]> result = new HashMap<Language,int[]>();
//    final int length = model.nSites();
//    int [] root = new int[length];
//    for (int i = 0; i < length; i++)
//      root[i] = SampleUtils.sampleMultinomial(rand, model.getInitialDistribution(i));
//    for (int c = 0; c < 2; c++)
//      _generate(model, depth-1, result, new Language(""+c), root, rand);
//    return result;
//  }
//  
//  private static void _generate(CTMC model, int depth,
//      Map<Language, int[]> result, Language curLang, int [] parent, Random rand)
//  {
//    // generate next
//    final int length = model.nSites();
//    int [] cur = new int[length];
//    for (int i = 0; i < length; i++)
//      cur[i] = SampleUtils.sampleMultinomial(rand, model.getTransitionPr(i,1.0)[parent[i]]);
//    if (depth==1)
//      result.put(curLang, cur);
//    else
//      for (int c = 0; c < 2; c++)
//        _generate(model, depth-1, result, new Language(curLang.toString()+c), cur, rand);
//  }
//  
//  private static CTMC generatingParams(int nSites)
//  {
//    double[][] Q = binaryRateMtx(1.0, new double[]{0.4,0.6});
//    return new CTMC.SimpleCTMC(Q, nSites);
//  }
//
////  public static Map<Language, String> lang2genus()//WalsData data)
////  {
////    Map<Language, String>  result = new HashMap<Language, String> ();
////    for (String lang : data.langWithMinNFeatures(0, null))
////      result.put(new Language(lang), data.getGenus(lang));
////    return result;
////  }
//  public static PartialCoalescentState initState(Map<Language, int[]> observations,
//      CTMC ctmc)
//  {
//    List<Language> leafNames = new ArrayList<Language>();
//    List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
//    for (Language lang : observations.keySet())
//    {
//      leafNames.add(lang);
//      leaves.add(DiscreteModelCalculator.observation(ctmc, observations.get(lang)));
//    }
//    return PartialCoalescentState.initialState(leaves, leafNames);
//  }
//  public static CTMC ctmc(double[][] statDistn, double rate)
//  {
//    List<double[][]> Qs = new ArrayList<double[][]>();
//    for (int site = 0; site < statDistn.length; site++)
//      Qs.add(binaryRateMtx(rate, statDistn[site]));
//    return new CTMC.GeneralCTMC(Qs);
//  }
//  public static double[][] binaryRateMtx(double rate, double[] stat)
//  {
//    if (!MathUtils.isProb(stat))
//      throw new RuntimeException();
//    double [][]result =  new double[][]{
//        {-rate*stat[1], rate*stat[1]},
//        { rate*stat[0],-rate*stat[0]}};
//    return result;
//  }
//  /**
//   * Laplace +1 estimator
//   * @param observations
//   * @return
//   */
//  private static double[][] statDistn(Map<Language, int[]> observations)
//  {
//    int nSites = observations.values().iterator().next().length;
//    double [][] result = new double[nSites][2];
//    for (int s = 0; s < nSites; s++)
//    {
//      for (Language lang : observations.keySet())
//      {
//        if (observations.get(lang)[s] != -1)
//          result[s][observations.get(lang)[s]]++;
//      }
//      result[s][0] += 1; result[s][1] += 1; // smooth
//      NumUtils.normalize(result[s]);
//    }
//    return result;
//  }
//
////  /**
////   * 1 - feature value is the most frequent
////   * 0 - otherwise
////   * @param data
////   * @param modes
////   * @param minNFeatures
////   * @param familyRestr
////   * @return
////   */
////  private static Map<Language, int[]> createObservations(WalsData data,
////      Map<Integer, Integer> modes, int minNFeatures, String familyRestr)
////  {
////    Map<Language, int[]> result = new HashMap<Language,int[]>();
////    List<Integer> features = new ArrayList<Integer>(data.allFeatures());
////    final int nSites = features.size();
////    for (String lang : data.langWithMinNFeatures(minNFeatures, familyRestr))
////    {
////      int [] current = new int[nSites];
////      for (int i = 0; i < nSites; i++)
////      {
////        Integer siteId = features.get(i);
////        if (!data.hasFeature(lang, siteId))
////          current[i] = -1;
////        else
////        {
////          Integer value = data.feature(lang, siteId);
////          if (value.equals(modes.get(siteId)))
////            current[i] = 1;
////          else
////            current[i] = 0;
////        }
////      }
////      result.put(new Language(lang), current);
////    }
////    return result;
////  }
//
////  /**
////   * 
////   * @param data
////   * @return site -> id of mode
////   */
////  private static Map<Integer, Integer> modes(WalsData data, int minNFeatures,String familyRestr)
////  {
////    Map<Integer,Integer> modes = new HashMap<Integer,Integer>();
////    for (Integer site : data.allFeatures())
////    {
////      Counter<Integer> counter = new Counter<Integer>();
////      for (String lang : data.langWithMinNFeatures(minNFeatures, familyRestr))
////        if (data.hasFeature(lang,site))
////          counter.incrementCount(data.feature(lang,site), 1.0);
////      modes.put(site, counter.argMax());
////    }
////    return modes;
////  }
//
//  public static void main(String [] args)
//  { 
//    ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
//    Execution.run(args, new WalsExperiment(pf), "wals", WalsDataset.class, "pf", pf);
//  }
}
