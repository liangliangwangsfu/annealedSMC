package pty.io;
import fig.basic.Option;
import fig.prob.SampleUtils;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import pty.learn.CTMCLoader;
import pty.smc.models.CTMC;
import pty.smc.models.CTMCUtils;

public class GeneratedDataset implements Dataset
{
  @Option public static int nSites = 100;
  @Option public static int depth = 4;
  /**
   * higher evaluationDepth means larger clusters
   */
  @Option public static int evaluationDepth = 2;
  @Option public static Random rand = new Random(1);

  public static CTMCLoader genCTMCLoader = new CTMCLoader();
  
  private CTMC ctmc = null;
  private Map<Taxon,int[]> generatedData = null;
  private Map<Taxon,String> populations = new HashMap<Taxon,String>();

  public Map<Taxon, String> getReferenceClusters()
  {
    return Collections.unmodifiableMap(populations);
  }

  public boolean hasReferenceClusters()
  {
    return true;
  }

  public int nCharacter(int site)
  {
    init();
    return ctmc.nCharacter(site);
  }

  public int nSites() { return nSites; }

  
  public Map<Taxon, double[][]> observations()
  {
    init();
    return Dataset.DatasetUtils.convert(Collections.unmodifiableMap(generatedData),this);
  }
  
  private Map<Taxon,int[]>  generate(CTMC model, int depth, Random rand)
  {
    if (depth <= 1) throw new RuntimeException();
    if (evaluationDepth > depth) throw new RuntimeException();
    Map<Taxon,int[]> result = new HashMap<Taxon,int[]>();
    final int length = model.nSites();
    int [] root = new int[length];
    for (int i = 0; i < length; i++)
      root[i] = SampleUtils.sampleMultinomial(rand, model.getInitialDistribution(i));
    for (int c = 0; c < 2; c++)
      _generate(model, depth-1, result, new Taxon(""+c), root, rand);
    return result;
  }
  private Set<Taxon> _generate(CTMC model, int depth,
      Map<Taxon, int[]> result, Taxon curLang, int [] parent, Random rand)
  {
    Set<Taxon> leaves = new HashSet<Taxon>();
    // generate next
    final int length = model.nSites();
    int [] cur = new int[length];
    for (int i = 0; i < length; i++)
      cur[i] = SampleUtils.sampleMultinomial(rand, model.getTransitionPr(i,1.0)[parent[i]]);  // TODO: fix this
    if (depth==1)
    {
      result.put(curLang, cur);
      leaves.add(curLang);
    }
    else
      for (int c = 0; c < 2; c++)
        leaves.addAll(_generate(model, depth-1, result, new Taxon(curLang.toString()+c), cur, rand));
    if (evaluationDepth == depth)
      for (Taxon lang : leaves)
        populations.put(lang, "Pop" + curLang.toString());
    return leaves;
  }
  
  private void init()
  {
    if (ctmc != null) return;
    genCTMCLoader.setData(new Dataset() {
      public Map<Taxon, String> getReferenceClusters()
      {
        throw new RuntimeException();
      }
      public boolean hasReferenceClusters()
      {
        return false;
      }
      public int nCharacter(int site)
      {
        throw new RuntimeException();
      }
      public int nSites()
      {
        return nSites;
      }
      public Map<Taxon, double[][]> observations()
      {
        throw new RuntimeException();
      }
    });
    ctmc = genCTMCLoader.load();
    CTMCUtils.saveInExec(ctmc, "generatingParams");
    generatedData = generate(ctmc, depth, rand);
  }
}
