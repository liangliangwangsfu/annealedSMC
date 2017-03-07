package pty.smc.test;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import ma.newick.NewickParser;
import nuts.util.Arbre;
import nuts.util.Tree;
import nuts.util.Arbre.ArbreMap;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import pty.RootedTree;
import pty.Train;
import pty.eval.Purity;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.LeaveOneOut;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.io.WalsDataset.WalsCorpusOperation;
import pty.learn.CTMCLoader;
import pty.learn.CTMCLoader.LoadingMethod;
import pty.smc.ConstrainedKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.PCSHash;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.models.BrownianModel;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.CTMC;
import pty.smc.test.SymmetryTest.MeanHeight;
import pty.smc.test.TestBrownianModel.KernelType;
import sage.LikelihoodModel;

public class MultiTests implements Runnable
{
  public void genTest(
      int nParticules, 
      int nSites, 
      int nRepeats, 
      KernelType kernelType, 
      Arbre<Taxon> trueTopo,
      double variance,
      int nThreads,
      DatasetType dataType,
      double trueH)
  {
    PartialCoalescentState init = null;
    PartialCoalescentState.disableUnrooted = true;
    Map<Taxon,String> refClusters = null;
    if (dataType == DatasetType.HGDP)
    {
//      if (siteSpec) throw new RuntimeException();
      // load the right amount of sites:
      HGDPDataset.maxNSites = nSites;
      HGDPDataset data = new HGDPDataset();
      final int actual = data.nSites();
      if (actual != nSites)
      { 
        LogInfo.warning("Requested " + nSites + " but there are only " + actual + " sites in data");
        nSites = actual;
      }
      BrownianModel model = new BrownianModel(data.nSites(), variance);
      init = null; //PartialCoalescentState.initState(data, model);
    }
    else if (dataType == DatasetType.WALS)
    {
      if (nSites!=Integer.MAX_VALUE) throw new RuntimeException();
      Dataset data = dataType.loadDataset();
      loader.rate = variance;
      loader.setData(data);
      CTMC ctmc = loader.load();
      init = PartialCoalescentState.initState(data, ctmc);
      refClusters = data.getReferenceClusters();
    }
    else throw new RuntimeException();
    // create the proper part. filter
    ParticleKernel<PartialCoalescentState> ppk = kernelType.load(init, null);
    ParticleFilter<PartialCoalescentState> pc = new ParticleFilter<PartialCoalescentState>();
    pc.N = nParticules;
    pc.nThreads = nThreads;
    pc.resampleLastRound = false;
    SummaryStatistics timeStats = new SummaryStatistics();
    SummaryStatistics perfStats = new SummaryStatistics();
    SummaryStatistics norm = new SummaryStatistics();
    SummaryStatistics hss = new SummaryStatistics();
    SummaryStatistics purity = new SummaryStatistics();
    SummaryStatistics loo = new SummaryStatistics();
    final int firstDatum = datum;
    for (int r = 0; r < nRepeats; r++)
    {
      ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> processor 
        = SymmetricDiff.createCladeProcessor();
      ParticleFilter.ParticleMapperProcessor<PartialCoalescentState,Set<Set<Taxon>>> mbr 
        = SymmetricDiff.createCladeProcessor();
      MeanHeight hpro
        = new MeanHeight();
      ParticleFilter.MAPDecoder<PartialCoalescentState> mapDecoder = new ParticleFilter.MAPDecoder<PartialCoalescentState>();
      ForkedProcessor<PartialCoalescentState> processors = new ForkedProcessor(processor, mapDecoder, mbr, hpro);
      // sample and time
      long time = System.currentTimeMillis();
      try {
        pc.sample(ppk,processors);
        time = System.currentTimeMillis() - time; 
        timeStats.addValue(time);
        // take various performance measurements
        Arbre<Taxon> recon=SymmetricDiff.clades2arbre(mbr.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE));
        RootedTree map = mapDecoder.map().getFullCoalescentState();
        if (dataType == DatasetType.WALS)
          loo.addValue(LeaveOneOut.loo(mapDecoder.map()));
        if (printTrees)
        {
          Train.outputTree(recon, "mbr-tree-datum-"+datum);
          Train.outputTree (map.topology(), "mapTree-" + datum, map.branchLengths());
        }
        if (trueTopo != null)
        {
          double perf = SymmetricDiff.normalizedSymmetricCladeDiff(recon, trueTopo);
          perfStats.addValue(perf);
          double reconH = hpro.ss.getMean();
          LogInfo.logs("!!!!!!!!!!!!!Recon=" + reconH);
          hss.addValue(Math.abs(trueH - reconH));
        }
        if (refClusters != null)
        {
          double pur = Purity.purity(Arbre.arbre2Tree(recon), refClusters);
          purity.addValue(pur);
        }
        norm.addValue(pc.estimateNormalizer());
        datum++;
      } catch (Exception e)
      {
        LogInfo.warning("One part. filter run failed");
        LogInfo.warning("Detail:" + e);
        e.printStackTrace();
        nRepeats++;
      }
    }
    LogInfo.logsForce(dataType+"-datum-" + firstDatum + "-" + (datum-1) + ": " + 
        "nPart="+nParticules+" "+
        "nSites="+nSites+" "+
        "nRep="+nRepeats+" "+
        "kernel="+kernelType+" "+
        "scale="+variance+" "+
        "nt="+nThreads+" "+
        "timeM="+timeStats.getMean()+" "+
        "timeSD="+timeStats.getStandardDeviation()+" "+
        "symmM="+(trueTopo==null?"N/A":perfStats.getMean())+" "+
        "symmSD="+(trueTopo==null?"N/A":perfStats.getStandardDeviation()) +" "+
        "hM="+(trueTopo==null?"N/A":hss.getMean())+" "+
        "hSD="+(trueTopo==null?"N/A":hss.getStandardDeviation()) +" "+
        "normM="+norm.getMean()+" "+
        "normSD=" + norm.getStandardDeviation()+ " "+
        "purM="+(purity==null?"N/A":purity.getMean()) + " " + 
        "purSD="+(purity==null?"N/A":purity.getStandardDeviation()) + " " +
        "looM="+loo.getMean() + " " +
        "looSD=" +loo.getStandardDeviation());

  }
  
  static int datum = 0;

  
  public void run () {
    if (nSitesIncr <= 0 || nSitesIncr <= 0)
      throw new RuntimeException();
    
    double trueH = Double.NaN;
    Arbre<Taxon> ar = null;
    try {
      NewickParser np = new NewickParser(IOUtils.openIn(trueTreePath));
      Tree<String> tree = np.parse();
      ar = Arbre.tree2Arbre(tree).preOrderMap(new ArbreMap<String,Taxon>() {
        @Override
        public Taxon map(Arbre<String> t) { return new Taxon(t.getContents()); }
      });
      trueH = RootedTree.Util.height(RootedTree.Util.load(new File(trueTreePath)));
      LogInfo.logs("True height:" + trueH);
    } catch (Exception e) { LogInfo.warning("Warning:"+trueTreePath+ " does not seem to be pointing to a " +
    		"valid topo...\nSkipping topo evaluation"); }
    
   
    for (int curScIt = 0; curScIt < scalingIters; curScIt++)
    {
      for (int curNSites = nSitesMin; curNSites <= nSitesMax && curNSites >=0; curNSites += nSitesIncr)
        for (int nParticles = nParticlesMin; nParticles <= nParticlesMax; nParticles+=nParticlesIncr)
          genTest(nParticles, curNSites, nRepeats, kernelType, ar, scaling, nThread,dataType,trueH); // just add purity here later
      scaling *= reScaling;
    }
    
    
  } 
  
  @Option public KernelType kernelType = KernelType.PRIOR_PRIOR;
  @Option public String trueTreePath = "none";
  @Option public int nRepeats = 10;
  @Option public int nSitesMin = Integer.MAX_VALUE;
  @Option public int nSitesMax = Integer.MAX_VALUE;
  @Option public int nSitesIncr = 100;
  @Option public int nThread = 8;
  @Option public double scaling = 1;
  @Option public int scalingIters = 1;
  @Option public double reScaling = 2;
  @Option public int nParticlesMin = 100;
  @Option public int nParticlesMax = 100;
  @Option public int nParticlesIncr = 100;
  @Option public DatasetType dataType = DatasetType.HGDP;
//  @Option public boolean siteSpec = false;
  
  @Option public boolean printTrees = false;
  
  public static CTMCLoader loader = new CTMCLoader();
  
  public static void main (String args []){
    if (Arrays.asList(args).contains("JARS"))
      Execution.jarFiles = new ArrayList<String>(Arrays.asList(
          "/home/eecs/bouchard/jars/ptychodus.jar",
          "/home/eecs/bouchard/jars/nuts.jar",
          "/home/eecs/bouchard/jars/fig.jar",
          "/home/eecs/bouchard/jars/pepper.jar"));
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    Execution.run(args, new MultiTests(),
        "hgdp", HGDPDataset.class,
        "wals", WalsDataset.class,
        "ctmc", loader,
        "bmc", BrownianModelCalculator.class);
  }
}
