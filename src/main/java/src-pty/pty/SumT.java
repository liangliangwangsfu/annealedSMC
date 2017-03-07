package pty;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import ma.newick.NewickParser;
import nuts.io.IO;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.EasyFormat;
import pty.eval.Purity;
import pty.io.LeaveOneOut;
import pty.io.WalsAnn;
import pty.io.WalsDataset;
import fig.basic.IOUtils;
import fig.basic.Interner;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.UnorderedPair;
import fig.exec.Execution;
import goblin.Taxon;


/**
 * WARNING: ASSUMES A ROOTED MODEL!!
 * @author bouchard
 *
 */
public class SumT  implements Runnable 
{
  @Option public String treeFile = null;
  @Option public boolean evaluatePurity = true;
  @Option public String purityGroupsFile = "";
  @Option public String fileFilter = "sample";
  @Option public int maxFileIndex = Integer.MAX_VALUE;
  @Option(gloss="fraction of files to be skipped after sorting and truncation") public double burnin = 0.1;
  @Option(gloss="number to skip between each tree processed") public int thinning = 0; 
  @Option public int evalFrequency = 0;
  
  private List<UnrootedTree> trees = new ArrayList<UnrootedTree>();
  private Map<Taxon,String> reference = null;
  private Counter<Set<Set<Taxon>>> topoCounter = new Counter<Set<Set<Taxon>>>();
  private Counter<Set<Taxon>> cladeCounter = new Counter<Set<Taxon>>();
  
  int loadedTreeCount = 0;
  public void loadTrees(File file)
  {
    LogInfo.track("Reading samples under " + file);
    List<File> subFiles =  new LinkedList<File>(IO.locate(file, IO.suffixFilter("gz","newick")));
    Iterator<File> iter = subFiles.iterator();
    while (iter.hasNext())
    {
      File cur = iter.next();
      if (!cur.getName().contains(fileFilter) ||
          indexExceedLimit(cur, maxFileIndex))
        iter.remove();
    }
    Collections.sort(subFiles, new Comparator<File>() {
      public int compare(File o1, File o2) {
        return new Integer(getIndex(o1)).compareTo(getIndex(o2));
      }
    });
    // compute the fraction of files to skip, and skip them
    int nFile = subFiles.size();
    int nToBurnIn = (int) (burnin * nFile);
    LogInfo.logsForce("Skipping " + nToBurnIn + " files (burned-in)");
    iter = subFiles.iterator();
    for (int i = 0; i < nToBurnIn; i++)
    {
      iter.next();
      iter.remove();
    }
    // process files
    LogInfo.logsForce("Processing files with thinning: " + thinning);
    int nTreeReadSoFar = 1;
    for (File subFile : subFiles)
    {
      BufferedReader reader = IOUtils.openInHard(subFile);
      String line = null;
      try {
        while ((line = reader.readLine()) != null)
          if (line.charAt(0) == '(')
          {
            if (thinning == 0 || (loadedTreeCount++)%thinning==0)
              trees.add(UnrootedTree.fromNewick(line));
          }
          else if (line.equals("")) ; 
          else
            LogInfo.warning("Warning skipped line in " + file.getAbsolutePath() + "\n\t" + line);
      reader.close();
      LogInfo.logs("" + (nTreeReadSoFar++) + "/" + subFiles.size() + ", " + trees.size() + " trees read so far");
      } catch (IOException e) 
      { LogInfo.warning("Problems with file:" + subFile); }
    }
    LogInfo.end_track();
  }

  private static boolean indexExceedLimit(File cur, int maxFileIndex)
  {
    return getIndex(cur) > maxFileIndex;
  }

  public static void main (String args [])
  {
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    Execution.run(args, new SumT(),
        "wals", WalsDataset.class);
  }
  
  private static class CladeInterner // only for memory purpose--dont want to risk with speed optimization
  {
    private Interner<Taxon> lInterner = new Interner<Taxon>();
    private Interner<Set<Taxon>> sInterner = new Interner<Set<Taxon>>();
    public Set<Set<Taxon>> internClades(Set<Set<Taxon>> clades)
    {
      Set<Set<Taxon>> result = new HashSet<Set<Taxon>>(clades.size());
      for (Set<Taxon> clade : clades)
        result.add(internClade(clade));
      return result;
    }
    private Set<Taxon> internClade(Set<Taxon> clade)
    {
      if (sInterner.isCanonical(clade))
        return clade;
      Set<Taxon> result = new HashSet<Taxon>(clade.size());
      for (Taxon lang : clade)
        result.add(lInterner.intern(lang));
      return sInterner.intern(result);
    }
  }

  @Override
  public void run() { consensus(); }
  public UnrootedTree consensus()
  {
    LogInfo.logs("");
    
    // load the wals reference if needed
    if (evaluatePurity)
    {
      if (purityGroupsFile.equals(""))
      {
        WalsDataset ds = WalsDataset.getPreprocessedCorpus();
        reference = ds.getReferenceClusters();
      }
      else
        reference = parseReferences(purityGroupsFile);
    }
    
    // load
//    for (String treeFile : treeFiles)
    loadTrees(new File(treeFile));
    LogInfo.logsForce("" + trees.size() + " trees read");
    // fill the counters
    LogInfo.track("Creating clade counters");
    CladeInterner interner = new CladeInterner();
    Counter<UnorderedPair<Taxon,Taxon>> pairwiseDistances =
      new Counter<UnorderedPair<Taxon,Taxon>>();
    SummaryStatistics hstat = new SummaryStatistics();
    int i = 1;
    for (UnrootedTree t : trees)
    {
      hstat.addValue(height(t));
      pairwiseDistances.incrementAll(t.pairwiseDistances());
      Set<Set<Taxon>> clades = t.clades();
      if (!topoCounter.containsKey(clades))
        clades = interner.internClades(clades);
      topoCounter.incrementCount(clades, 1.0);
      for (Set<Taxon> clade : clades)
        cladeCounter.incrementCount(clade, 1.0);
      LogInfo.logs("" + (i++) + "/" + trees.size() + " [" + cladeCounter.size() + " unique clades]");
      if (evaluatePurity && evalFrequency != 0 && (i % evalFrequency == 0))
      {
        LogInfo.logsForce("Purity of current sample (" + i + "): " + Purity.purity(t.clades(), reference));
        LogInfo.logsForce("MBR Purity so far: " + Purity.purity(findMin(topoCounter, cladeCounter), reference));
      }
    }
    LogInfo.end_track();
    
    LogInfo.logs("Approx reconstructed dist:" + hstat.getMean());

    double norm = trees.size();
    for (UnorderedPair<Taxon,Taxon> key : pairwiseDistances.keySet())
      pairwiseDistances.setCount(key, pairwiseDistances.getCount(key) / norm);
    
    final String distanceFile = Execution.getFile("distances"); 
    IO.writeToDisk(distanceFile, phylipDistanceMatrix(pairwiseDistances));

    Set<Set<Taxon>> recontructedClades = null;
    // 3: reconstruct
    if (trees.size() == 1)
      recontructedClades = trees.get(0).clades();
    else
    {
      LogInfo.track("MBR reconstruction", false);
      LogInfo.logs("Number of unique topologies: " + topoCounter.size());
      recontructedClades = findMin(topoCounter, cladeCounter);
      if (trees.size() != (int) topoCounter.totalCount())
        throw new RuntimeException(trees.size() + " vs " + (int) topoCounter.totalCount());
      LogInfo.end_track();
    }
    LogInfo.logsForce("MBR clades: " + recontructedClades);
    // output ranked clades
    int rk = 1;
    LogInfo.track("Clades:", true);
    for (Set<Taxon> clade : cladeCounter)
      if (recontructedClades.contains(clade))
        if (clade.size() > 1 && clade.size() < nLeaves(recontructedClades) - 1)
          LogInfo.logs("" + (rk++) + "\t" + (cladeCounter.getCount(clade)/trees.size()) + "\t" +clade);
    LogInfo.end_track();
    // 4: write purity
    if (evaluatePurity)
    {
      LogInfo.track("Evaluation",true);
//      LogInfo.logs("Reference: " + reference);
      LogInfo.logs("Purity: " + Purity.purity(recontructedClades, reference));
      LogInfo.end_track();
    }
    // output an argmax tree
    UnrootedTree result = null;
    for (UnrootedTree nct : trees)
      if (nct.clades().equals(recontructedClades))
      {
        result = nct;
        IO.writeToDisk(new File(Execution.getFile("consensus.newick")), nct.toNewick());
        break;
      }
    // reconstruct using pairwise distances and phylip
    if (evaluatePurity)
    {
      String recon = IO.call("/bin/bash cmds/phylipnj.bash " + distanceFile);
      UnrootedTree pairRecon = UnrootedTree.fromNewick(recon);
      LogInfo.logs("Purity (from distances): " + Purity.purity(pairRecon.clades(), reference));
    }
    return result;
  }
  private double height(UnrootedTree t)
  {
    double min = Double.POSITIVE_INFINITY;
    // find the center
    for (Taxon lang : t.getTopology().vertexSet())
      if (!t.leavesSet().contains(lang))
      {
        double cur = meanD(lang, t);
        if (cur < min)
          min = cur;
      }
    return min;
  }

  private double meanD(Taxon lang, UnrootedTree t)
  {
    double mean = 0.0;
    for (Taxon leaf : t.leaves())
      mean += t.totalBranchLengthDistance(leaf, lang);
    return mean / ((double) t.leaves().size());
  }

  public static String phylipDistanceMatrix(
      Counter<UnorderedPair<Taxon, Taxon>> pairwiseDistances)
  {
    StringBuilder result = new StringBuilder();
    Set<Taxon> langsSet = CollUtils.set();
    for (UnorderedPair<Taxon, Taxon> key : pairwiseDistances.keySet())
    {
      langsSet.add(key.getFirst()); langsSet.add(key.getSecond());
    }
    List<Taxon> langs = new ArrayList<Taxon>(langsSet);
    result.append(langs.size() + "\n");
    for (Taxon lang : langs)
    {
      result.append(WalsAnn.cleanForPhylip(lang.toString()));
      for (Taxon lang2 : langs)
        result.append("  " + EasyFormat.fmt(pairwiseDistances.getCount(new UnorderedPair<Taxon,Taxon>(lang,lang2))));
      result.append('\n');
    }
    return result.toString();
  }

  public static Map<Taxon, String> parseReferences(String purityGroupsFile2)
  {
    Map<Taxon, String> result = CollUtils.map();
    for (String line : IO.i(purityGroupsFile2))
    {
      String [] fields = line.split("\\s+");
      result.put(new Taxon(fields[0]), fields[1]);
    }
    return result;
  }

  private int nLeaves(Set<Set<Taxon>> clades)
  {
    int max = 0; 
    for (Set<Taxon> clade : clades)
      if (clade.size() > max)
        max = clade.size();
    return max;
  }

  public static Set<Set<Taxon>> findMin(
      Counter<Set<Set<Taxon>>> topoCounter,
      Counter<Set<Taxon>> cladeCounter)
  {
    Set<Set<Taxon>> argmin = null;
    double minValue = Double.POSITIVE_INFINITY;
    int i = 1;
    for (Set<Set<Taxon>> currentElt : topoCounter.keySet())
    {
      // old method
//      int currentValue = sumOfLosses(currentElt, topoCounter, SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE);
      final int currentValue = sumOfSymmLosses(currentElt, cladeCounter, (int)topoCounter.totalCount());
//      if (currentValue != newMethod)
//        System.out.println("CV=" + currentValue + ",NM=" + newMethod);
      if (currentValue < minValue)
      {
        minValue = currentValue;
        argmin = currentElt;
      }
//      LogInfo.logs("" + (i++) + "/" +  topoCounter.keySet().size());
    }
    return argmin;
  }
  
  private static int sumOfSymmLosses(Set<Set<Taxon>> currentElt, Counter<Set<Taxon>> cladeCounts, int nTrees)
  {
    if (nTrees == 1) throw new RuntimeException();
    final int  nCladesPerTree = currentElt.size();
//    System.out.println(nCladesPerTree);
    int result = (nTrees - 1) * nCladesPerTree;
    for (Set<Taxon> clade : currentElt)
      result = result - ((int)cladeCounts.getCount(clade) - 1);
    return 2*result;
  }
  
  public static int getIndex(File f)
  {
    try { return Integer.parseInt(f.getName().replaceAll("[^0-9]*","")); }
    catch (Exception e) { return -1; }
  }

//  private static <T> int sumOfLosses(T currentElt, Counter<T> elementMultiplicities, LossFct<T> lossFct)
//  {
//    int sum = 0;
//    for (T otherElement : elementMultiplicities.keySet())
//      sum += elementMultiplicities.getCount(otherElement) * lossFct.loss(currentElt, otherElement);
//    return sum;
//  }

}
