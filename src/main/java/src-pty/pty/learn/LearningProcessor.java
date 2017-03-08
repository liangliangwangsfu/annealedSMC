/**
 * 
 */
package pty.learn;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;

import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.Parallelizer;
import goblin.Taxon;
import nuts.math.GMFct;
import nuts.math.RateMtxUtils;
import nuts.tui.Table;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import pty.RootedTree;
import pty.Observations;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.models.CTMC;
import Jama.Matrix;

public class LearningProcessor implements ParticleProcessor<PartialCoalescentState>
{
  @Option public static int nThreads = 1;
  @Option public static boolean useReversibleMtx = false;
  @Option public static boolean tieSites = true;
  
  public static Matrix getSufficientStatistics(
      RootedTree state,
      CTMC ctmc,
      Observations observations,
      int site)
  {
    final GMFct<Taxon> posterior = DiscreteBP.posteriorMarginalTransitions(state, ctmc, observations, site);
    final int nCharacters = ctmc.nCharacter(site);
    final Matrix result = new Matrix(nCharacters, nCharacters);
    for (Arbre<Taxon> node : state.topology().nodes())
      if (!node.isRoot())
      {
        Taxon par = node.getParent().getContents(), 
        cur = node.getContents();
        double T = state.branchLengths().get(cur);
        double [][][][] expectations = CTMCExpectations.expectations(T, ctmc.getRateMtx(site));
        for (int topState = 0; topState < nCharacters; topState++)
         for (int botState = 0; botState < nCharacters; botState++)
         {
           Matrix current = new Matrix(expectations[topState][botState]);
           double post = posterior.get(par, cur, topState, botState);
           result.plusEquals(current.times(post));
         }
      }
    return result;
  }
//  public static Matrix getSufficientStatistics(PartialCoalescentState pcs, final int site)
//  {
//    Pair<Arbre<Language>,Map<Language,Double>> p = pcs.getArbreAndBranchLengths();
//    final CTMC ctmc = pcs.getCTMC();
//    Map<Language,double[]> obs = pcs.getObservations(site);
//    return getSufficientStatistics(p.getFirst(), p.getSecond(), ctmc, obs, site);
//    final GMFct<Language> posterior = DiscreteBP.posteriorMarginalTransitions(pcs,site);
//    Pair<Arbre<Language>,Map<Language,Double>> p = pcs.getArbreAndBranchLengths();
//    final CTMC ctmc = pcs.getCTMC();
//    final int nCharacters = ctmc.nCharacter(site);
//    final Matrix result = new Matrix(nCharacters, nCharacters);
//    Arbre<Language> a = p.getFirst();
//    final Map<Language,Double> bls = p.getSecond();
//    for (Arbre<Language> node : a.nodes())
//      if (!node.isRoot())
//      {
//        Language par = node.getParent().getContents(), 
//        cur = node.getContents();
//        double T = bls.get(cur);
//        double [][][][] expectations = CTMCExpectations.expectations(T, ctmc.getRateMtx(site));
//        for (int topState = 0; topState < nCharacters; topState++)
//         for (int botState = 0; botState < nCharacters; botState++)
//         {
//           Matrix current = new Matrix(expectations[topState][botState]);
//           double post = posterior.get(par, cur, topState, botState);
//           result.plusEquals(current.times(post));
//         }
//      }
//    return result;
//  }
  
  private final Matrix [] suffStats;  // one for each site
  public LearningProcessor(CTMC current) 
  {
    this.suffStats = new Matrix[current.nSites()];
    for (int s =0 ; s < current.nSites();s++)
      suffStats[s] = new Matrix(current.nCharacter(s),current.nCharacter(s));
  }
  public void process(final PartialCoalescentState pcs, final double weight)
  {
    process(pcs.getFullCoalescentState(), pcs.getCTMC(), pcs.getObservations(), weight);
//    Pair<Arbre<Language>,Map<Language,Double>> p = pcs.getArbreAndBranchLengths();
//    final CTMC ctmc = pcs.getCTMC();
//    Map<Language,double[]> obs = pcs.getObservations(site);
//    return getSufficientStatistics(p.getFirst(), p.getSecond(), ctmc, obs, site);
  }
  public void process(      
      final RootedTree state,
      final CTMC ctmc,
      final Observations observations,
      final double weight)
  {
    List<Integer> ints = CollUtils.ints(suffStats.length);
    Parallelizer<Integer> par = new Parallelizer<Integer>(nThreads);
    par.setPrimaryThread();
    par.process(ints, new Parallelizer.Processor<Integer>() {
      public void process(Integer s, int i, int n, boolean log)
      {
        Matrix current = getSufficientStatistics(state,ctmc, observations, s).times(weight);
        synchronized (suffStats) { suffStats[s].plusEquals(current); }
      }
    });
  }
  public CTMC reestimate(CTMC old)
  {
    if (useReversibleMtx)
    {
      if (!tieSites) 
        throw new RuntimeException();
      if (!old.isSiteTied()) 
        throw new RuntimeException();
      double [] sd = old.getInitialDistribution(0);
      return getTiedReversibleMLE(sd);
    }
    else 
    {
      if (tieSites)
        return getTiedMLE();
      else 
        return getMLE();
    }
  }
  private CTMC getMLE()
  {
    List<double[][]> Qs = new ArrayList<double[][]>();
    for (Matrix m :  suffStats)
      Qs.add(Estimators.getGeneralRateMatrixMLE(m));
    return new CTMC.GeneralCTMC(Qs);
  }
  public CTMC getTiedMLE()
  {
    return new CTMC.SimpleCTMC(Estimators.getGeneralRateMatrixMLE(sumSuffStats()),suffStats.length);
  }
  public Matrix sumSuffStats()
  {
    int nCharacters = suffStats[0].getColumnDimension();
    Matrix sum = new Matrix(nCharacters, nCharacters);
    for (Matrix m : suffStats)
    {
      if (m.getColumnDimension() != nCharacters)
        throw new RuntimeException();
      sum.plusEquals(m);
    }
    return sum;
  }
  /**
   * assumes uniform stat. distn
   * @return
   */
  public CTMC getTiedReversibleMLE()
  {
    int nChars = suffStats[0].getColumnDimension();
    double [] statDistn = new double[nChars];
    for (int i = 0; i < statDistn.length; i++)
      statDistn[i] = 1.0 / ((double)nChars);
    return getTiedReversibleMLE(statDistn);
  }
  public CTMC getTiedReversibleMLE(double [] statDistn)
  {
    int nCharacters = suffStats[0].getColumnDimension();
    double [][] rateMtx = new double[nCharacters][nCharacters];
    Matrix sum = sumSuffStats();
//    try{
    for (int r = 0; r < nCharacters; r++)
      for (int c = 0; c < nCharacters; c++)
        if (r != c)
          rateMtx[r][c] = statDistn[c] * (sum.get(r,c)+sum.get(c,r)) 
            / (sum.get(c,c)*statDistn[r] + sum.get(r,r)*statDistn[c]);
    RateMtxUtils.fillRateMatrixDiagonalEntries(rateMtx);
    return new CTMC.SimpleCTMC(rateMtx, suffStats.length);
//    }
//    catch (Exception e)
//    {
//      e.printStackTrace();
//      LogInfo.logs("stat dist:" + Arrays.toString(statDistn)); //////////////////
//      LogInfo.logs("sum suff stats:\n" + Table.toString(sum)); 
//      throw new RuntimeException();
//    }
  }
  @Override
  public String toString()
  {
    StringBuilder result = new StringBuilder();
    for (Matrix m : suffStats)
      result.append(Table.toString(m)+'\n');
    return result.toString();
  }
}