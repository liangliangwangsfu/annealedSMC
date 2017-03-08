package pty.learn;
import java.io.*;
import java.util.*;

import nuts.lang.ArrayUtils;
import nuts.math.RateMtxUtils;
import nuts.tui.Table;
import nuts.util.CollUtils;
import nuts.util.Indexer;
import nuts.util.MathUtils;

import pty.io.Dataset;
import pty.io.WalsDataset;
import pty.io.WalsProcessingScript;
import pty.io.WalsDataset.BioCharacter;
import pty.io.WalsDataset.Site;
import pty.smc.models.CTMC;
import pty.smc.models.CTMCUtils;

import ma.SequenceType;

import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import goblin.Taxon;

/**
 * Different ways of loading a CTMC model
 * @author bouchard
 *
 */
public class CTMCLoader
{
  public static enum LoadingMethod 
  {
    BUILT_IN {
      @Override
      public CTMC load(CTMCLoader loader)
      {
        if (loader.siteSpecific) 
          throw new RuntimeException();
        CTMC result = CTMC.SimpleCTMC.fromSequenceType(loader.data.nSites(), loader.builtInSequenceType, loader.rate);
        LogInfo.logs("Simple built-in CTMC:\n" + result);
        return result;
      }
    },
    SCRIPT {
      @Override
      public CTMC load(CTMCLoader loader)
      {
        if (loader.siteSpecific) 
          throw new RuntimeException();
        WalsDataset dataset = (WalsDataset) (loader.data);
        WalsProcessingScript script = WalsDataset.getScript();
        List<double[][]> Qs = CollUtils.list();
        for (int s = 0; s < dataset.nSites(); s++)
        {
          final Site site = dataset.siteIndexer().i2o(s);
          final int nChars = dataset.nCharacter(s);
          double [][] current = new double[nChars][nChars];
          // need to make sure the characters are indexed in increasing order
          Indexer<BioCharacter> characterIndexer = dataset.charIndexers().get(site);
          // determine if this site is ordered & create the matrix
          if (script.orderFeature.contains(site.toString()))
          {
            // create the matrix
            for (int c1 = 0 ; c1 < nChars; c1++)
              // create matrix
              for (int c2 = 0; c2 < nChars; c2++)
                if (c1 != c2 && linked(characterIndexer, c1, c2))
//                if (Math.abs(characterIndexer.i2o(c1).index-characterIndexer.i2o(c2).index) == 1)
                  current[c1][c2] = 1.0;
          }
          else
            for (int c1 = 0; c1 < nChars; c1++)
              for (int c2 = 0; c2 < nChars ; c2++)
                if (c1 != c2)
                  current[c1][c2] = 1.0;
          RateMtxUtils.fillRateMatrixDiagonalEntries(current);
          Qs.add(current);
        }
        return new CTMC.GeneralCTMC(Qs);
      }

      
    },
    FILE {
      @Override
      public CTMC load(CTMCLoader loader)
      {
        return CTMCUtils.unSerialize(new File(loader.file));
      }
    },
    HEURISTIC_ESTIMATE {
      @Override
      public CTMC load(CTMCLoader loader)
      {
        return loader.getHeuristicEstimate();
      }
    };
    public abstract CTMC load(CTMCLoader loader);
  }
  @Option public LoadingMethod loadingMethod = LoadingMethod.HEURISTIC_ESTIMATE;
  
  // specific to FILE
  @Option public String file = "init.CTMC";
  
  // specific to BUILT_IN
  @Option public SequenceType builtInSequenceType = SequenceType.BINARY;
  
  // for BUILT_IN and HEURISTIC_ESTIMATE
  @Option public double rate = 1.0;
  
  // specific to HEURISTIC_ESTIMATE, 
  @Option public boolean siteSpecific = false;
  @Option public boolean forceUniform = false;
  
  public CTMC load() { return loadingMethod.load(this); }
  public Dataset data;
  public void setData(Dataset data) { this.data = data; }
  
  public static boolean linked(Indexer<BioCharacter> characterIndexer, int c1,
      int c2)
  {
//    if (characterIndexer.size() == 6)
//      System.out.println("fuck");
    // check that c2 is that c1 and c2 are closest in the original indexing
    final int _c1Index = characterIndexer.i2o(c1).index;
    final int _c2Index = characterIndexer.i2o(c2).index;
    // wlog, c1 < c2
    final int c1Index = _c1Index < _c2Index ? _c1Index : _c2Index;
    final int c2Index = _c1Index < _c2Index ? _c2Index : _c1Index;
    for (int c3 = 0; c3 < characterIndexer.size(); c3++)
    {
      final int c3Index = characterIndexer.i2o(c3).index;
      if (c1Index < c3Index && c3Index < c2Index)
        return false;
    }
    return true;
  }
  
  private CTMC getHeuristicEstimate()
  {
    if (siteSpecific)
    {
      List<double[][]> rateMtrices = new ArrayList<double[][]>();
      double [][] sds = statDistn(data.observations());
      for (int s = 0; s < sds.length; s++)
        rateMtrices.add(RateMtxUtils.reversibleRateMtx(rate, sds[s]));
      return new CTMC.GeneralCTMC(rateMtrices);
    }
    else
    {
      double [] sd = globalStatDistn(data.observations());
      double [][] rateMtx = RateMtxUtils.reversibleRateMtx(rate, sd);
      LogInfo.logs("Estimated stat dist:" + Arrays.toString(sd));
      LogInfo.logs("Estimated rate matrix from stat dist:\n" + Table.toString(rateMtx));
      return new CTMC.SimpleCTMC(rateMtx, data.nSites());
    }
  }
  
  /**
   * Laplace +1 estimator
   * @param observations
   * @return
   */
  private double[][] statDistn(Map<Taxon, double[][]> observations)
  {
    int nSites = data.nSites();
    double [][] result = new double[data.nSites()][];
    if (forceUniform)
    {
      double [][] anObsArray = observations.values().iterator().next();
      for (int s = 0; s < nSites; s++)
      {
        int nChars = anObsArray[s].length;
        result[s] = new double[nChars];
        for (int i = 0; i < result[s].length; i++)
          result[s][i] = 1.0/((double)nChars);
      }
      return result;
    }
    for (int s = 0; s < nSites; s++)
    {
      result[s] = new double[data.nCharacter(s)];
      for (Taxon lang : observations.keySet())
//        if (observations.get(lang)[s] != Dataset.unknownCode)
          process(observations.get(lang)[s], result[s]);
      for (int c = 0; c < data.nCharacter(s); c++)
        result[s][c] += 1;  // smooth
      NumUtils.normalize(result[s]);
    }
    return result;
  }
  private void process(double currentObs[], double [] result)
  {
    double sum = MathUtils.sum(currentObs);
    if (sum == currentObs.length) return;
    else if (MathUtils.close(1.0, sum))
      for (int i =0 ; i < currentObs.length; i++)
        result[i] += currentObs[i];
    else
      throw new RuntimeException();
//    if (currentObs > 0 && currentObs < 1)
//    {
//      result[0] += currentObs;
//      result[1] += 1.0 - currentObs;
//    }
//    else if (((int) currentObs) - currentObs == 0)
//      result[(int)currentObs]++;  // TODO: fix this
//    else
//      throw new RuntimeException();
  }
  private double[] globalStatDistn(Map<Taxon, double[][]> observations)
  {
    if (forceUniform) throw new RuntimeException();
    int nChars = data.nCharacter(0);
    double [] result = new double[nChars];
    for (Taxon lang : observations.keySet())
      for (int s = 0; s < data.nSites(); s++)
//        if (observations.get(lang)[s] != Dataset.unknownCode)
          process(observations.get(lang)[s], result);
    NumUtils.normalize(result);
    return result;
  }
}
