package pty.smc.exp;
import java.io.*;
import java.util.*;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.SumT;
import pty.eval.SymmetricDiff;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.mcmc.Main;
import pty.mcmc.ParallelTemperingChain;
import pty.mcmc.PhyloSampler;
import pty.mcmc.ProposalDistribution;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.NoLikelihoodModel;
import pty.smc.test.GenerateFromTree;
import pty.smc.test.MultiTests;

import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.DataPrepUtils;
import goblin.Taxon;

import nuts.io.IO;
import nuts.util.Arbre;
import nuts.util.Tree;

public class TimingExp implements Runnable
{
  private static GenerateFromTree genFromTree = new GenerateFromTree();
  private static MultiTests pfTest = new MultiTests();
  private static Main mcmcMain = new Main();
  
  @Option public boolean forceVarianceEqualToGenerated = true;
  @Option public Random mcmcInitRandom = new Random(1);
  @Option public double ratioOfTimeForMCMC = 1.0; // set this to take decoding into account
  @Option public boolean runSMC = true;
  @Option public double timeIfMCMCOnly = 10000;

  /**
   * @param args
   */
  public static void main(String[] args)
  {
    IO.run(args, new TimingExp(), 
        "mt", pfTest,
        "gft", genFromTree,
        "bmc", BrownianModelCalculator.class,
        "mcmcMain", mcmcMain, 
        "prior", PhyloSampler._defaultPriorOptions,
        "prop", ProposalDistribution.Util._defaultProposalDistributionOptions,
        "sampler", PhyloSampler._defaultPhyloSamplerOptions,
        "partemp", ParallelTemperingChain._defaultTemperingOptions);
  }
  
  public static Tree<Taxon> removeUselessInternalNodes(Tree<Taxon> t)
  {
    List<Tree<Taxon>> newChildren = new ArrayList<Tree<Taxon>>();
    if (t.getChildren().size() == 1) return removeUselessInternalNodes(t.getChildren().get(0));
    else
      for (Tree<Taxon> child : t.getChildren())
        newChildren.add(removeUselessInternalNodes(child));
    return new Tree<Taxon>(t.getLabel(), newChildren);
  }

  @Override
  public void run()
  {
    try{
      RootedTree goldTree = RootedTree.Util.load(new File(genFromTree.treeFile));
      // generate some data according to this tree
      LogInfo.track("Generating data");
      File data = genFromTree.generate();
      LogInfo.end_track();
      
      // FORCE VARIANCES
      LogInfo.logsForce("Forcing variance of the pf test to be equal to the generated one: " + genFromTree.scale);
      pfTest.scaling = genFromTree.scale;
      pfTest.scalingIters = 1;
      pfTest.dataType = DatasetType.HGDP;
      
      long time = (long) timeIfMCMCOnly;
      
      if (runSMC)
      {
        // run the particle filter on that data several times
        pfTest.trueTreePath = genFromTree.treeFile;
        HGDPDataset.path = data.getAbsolutePath();
        LogInfo.track("Running pf");
        time = System.currentTimeMillis();
        pfTest.run();
        time = (System.currentTimeMillis() - time) / pfTest.nRepeats;
        LogInfo.logsForce("Mean time for PF: " + time);
        LogInfo.end_track();
      }

      LogInfo.track("Running mcmc");
      
      SummaryStatistics diffs = new SummaryStatistics();
      SummaryStatistics timeStat = new SummaryStatistics();
      for (int curRep = 0; curRep < pfTest.nRepeats; curRep++)
      {
        long mcmcTime = System.currentTimeMillis();
        // create a random init for MCMC
        UnrootedTree mcmcInit = randomInit(UnrootedTree.fromNewick(new File(genFromTree.treeFile)).leaves(), mcmcInitRandom);
        File mcmcInitFile = new File(Execution.getFile("mcmcInit-" + curRep + ".newick"));
        IO.writeToDisk(mcmcInitFile, mcmcInit.toNewick());
        
        // run the MCMC for as much time as pf took in average
        ParallelTemperingChain._defaultTemperingOptions.timeCutOff = (long) (time * ratioOfTimeForMCMC);
        mcmcMain.dataTypes = new ArrayList<DatasetType>(Arrays.asList(DatasetType.HGDP));
        HGDPDataset.path = data.getAbsolutePath();
        mcmcMain.brownianMotionVariance = genFromTree.scale;
        mcmcMain.pathsToInitTree = new ArrayList<String>(Arrays.asList(mcmcInitFile.getAbsolutePath()));
        mcmcMain.conditionOnFamilies = false;
        PhyloSampler._defaultPhyloSamplerOptions.logFrequency = 1;
        
        LogInfo.track("MCMC run " + curRep + "/" + pfTest.nRepeats);
        File samples = mcmcMain.run("repeat" + curRep).get(0);
        LogInfo.end_track();
        
        // decode
        SumT sumT = new SumT();
        sumT.treeFile = samples.getAbsolutePath();
        sumT.evaluatePurity = false;
        UnrootedTree consensus = sumT.consensus();
        
        // evaluate by picking the oracle rooting
        diffs.addValue(eval(consensus, goldTree));
        timeStat.addValue(System.currentTimeMillis() - mcmcTime);
      }
      // print average evaluation
      LogInfo.logsForce("Mean time for MCMC: " + timeStat.getMean());
      LogInfo.logsForce("SD time for MCMC: " + timeStat.getStandardDeviation());
      LogInfo.logsForce("Mean symm diff for MCMC: " + diffs.getMean());
      LogInfo.logsForce("SD of symm diff for MCMC: " + diffs.getStandardDeviation());
      LogInfo.end_track();
      
    }catch (Exception e) { throw new RuntimeException(e); }
  }

  private double eval(UnrootedTree consensus, RootedTree goldTree)
  {
    Arbre<Taxon> ref = goldTree.topology();
    double min = Double.POSITIVE_INFINITY;
    for (Tree<Taxon> lang : consensus.toTree().getPostOrderTraversal())
    {
      Tree<Taxon> rooted = consensus.toTree(lang.getLabel());
      double cur = SymmetricDiff.normalizedSymmetricCladeDiff(Arbre.tree2Arbre(rooted), ref);
      if (cur < min) min = cur; 
    }
    return min;
  }

  private UnrootedTree randomInit(List<Taxon> leaves, Random rand)
  {
    LinkedList<Tree<Taxon>> subtrees = new LinkedList<Tree<Taxon>>();
    for (Taxon leaf : leaves)
      subtrees.add(new Tree<Taxon>(leaf));
    int id = 0;
    while (subtrees.size() > 3)
    {
      Tree<Taxon> 
      subt1 = subtrees.remove(rand.nextInt(subtrees.size())),
      subt2 = subtrees.remove(rand.nextInt(subtrees.size())); 
      Tree<Taxon> newT = new Tree<Taxon>(new Taxon("intern_" + (id++)));
      newT.getChildren().add(subt1);
      newT.getChildren().add(subt2);
      subtrees.add(newT);
    }
    if (subtrees.size() != 3) 
      throw new RuntimeException();
    Tree<Taxon> resultTopo = new Tree<Taxon>(new Taxon("intern_" + (id++)));
    for (Tree<Taxon> t : subtrees)
      resultTopo.getChildren().add(t);
    Map<Taxon,Double> bl = new HashMap<Taxon,Double>();
    for (Tree<Taxon> lang : resultTopo.getPostOrderTraversal())
      bl.put(lang.getLabel(), rand.nextDouble());
    return UnrootedTree.fromNewick(DataPrepUtils.newick(resultTopo, bl, true));
  }

}
