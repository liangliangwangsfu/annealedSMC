package pty;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ma.newick.NewickParser;
import nuts.io.IO;
import nuts.util.Arbre;
import nuts.util.Tree;
import nuts.util.Arbre.ArbreMap;
import pty.eval.Purity;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.GeneratedDataset;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.learn.CTMCLoader;
import pty.learn.DiscreteBP;
import pty.learn.LearningProcessor;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPostKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.CTMCUtils;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.DataPrepUtils;
import goblin.Taxon;

public class Train implements Runnable
{
  @Option public int nEMIters = 20;
  @Option public DatasetType datasetType = DatasetType.WALS;

  @Option public int increaseNSamplesPerEMIteration = 0;
  @Option public EstimationMethod estimationMethod = EstimationMethod.UNSUPERVISED;
  
  @Option public String guideTreePath = autoGuideTree;
  public static final String autoGuideTree = "AUTO";
  @Option public double guideTreeScaling = 100;
  
  @Option public boolean usePriorPost = false;
  
//  @Option public boolean scoreAllSamples;
  
  public static enum EstimationMethod { NONE, UNSUPERVISED, SUPERVISED; }
  
  private static ParticleFilter<PartialCoalescentState> pf;
  private static CTMCLoader loader;
  
  private Dataset data;
  
  public void run()
  {
    // load data, init params
    data = datasetType.loadDataset(); 
    loader.setData(data);
    CTMC ctmc = loader.load();
    CTMCUtils.saveInExec(ctmc, "init");
    // do EM
    for (int emIter = 0; emIter < nEMIters; emIter++)
    {
      LogInfo.track("EM Iteration " + (emIter+1) + "/" + nEMIters,true);
      // prepare sampling machinery
      PartialCoalescentState initState = PartialCoalescentState.initState(data, ctmc);
      ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> mbr =
        SymmetricDiff.createCladeProcessor();
      LearningProcessor ssp = new LearningProcessor(ctmc);
//      ParticleMapperProcessor<PartialCoalescentState,PartialCoalescentState> pmp =
//        ParticleMapperProcessor.saveParticlesProcessor();
      ForkedProcessor<PartialCoalescentState> processors 
        = new ForkedProcessor<PartialCoalescentState>(mbr);
//      if (scoreAllSamples && data.hasReferenceClusters())
//        processors.processors.add(pmp);
      if (estimationMethod == EstimationMethod.UNSUPERVISED) 
        processors.processors.add(ssp);
      ParticleKernel<PartialCoalescentState> kernel = 
        (usePriorPost ? new PriorPostKernel(initState) :  new PriorPriorKernel(initState));
      // E step : sample!
      pf.sample(kernel, processors);
      pf.N += increaseNSamplesPerEMIteration;
      // reconstruct a tree
      Arbre<Taxon> reconstruction = outputTree( SymmetricDiff.clades2arbre(mbr.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE)), "consensusTree-" + emIter);
      // evaluate purity if possible
      if (data.hasReferenceClusters())
      {
        LogInfo.logs("Purity:" + 
            Purity.purity(Arbre.arbre2Tree(reconstruction), 
                          data.getReferenceClusters()));
//        if (scoreAllSamples)
//          scoreAllSamples(pmp.getCounter(), data.getReferenceClusters(),
//              reconstruction);
      }
      // M step : update parameters
      if (estimationMethod == EstimationMethod.UNSUPERVISED)
      {
        LogInfo.track("Unsupervised reestimation of parameters");
        ctmc = ssp.reestimate(ctmc);
        CTMCUtils.saveInExec(ctmc, "unsup-reest-" + (emIter+1));
        LogInfo.end_track();
      }
      else if (estimationMethod == EstimationMethod.SUPERVISED)
      {
        LogInfo.track("Supervised reestimation of parameters");
        // a single, guide tree is processed
        ssp.process(getGuideCoalescent(guideTreePath,data,guideTreeScaling), ctmc, data, 1.0); 
        LogInfo.logs("Data likelihood before reestimation: " + 
            DiscreteBP.dataLogLikelihood(getGuideCoalescent(guideTreePath,data,guideTreeScaling), ctmc, data));
        ctmc = ssp.reestimate(ctmc);
        LogInfo.logs("Data likelihood after reestimation: " + 
            DiscreteBP.dataLogLikelihood(getGuideCoalescent(guideTreePath,data,guideTreeScaling), ctmc, data));
        CTMCUtils.saveInExec(ctmc, "sup-reest-" + (emIter+1));
        LogInfo.end_track();
      }
      LogInfo.end_track();
    }
  }
  
//  private void scoreAllSamples(
//      Counter<PartialCoalescentState> particles,
//      Map<Language, String> referenceClusters,
//      Arbre<Language> _reconstruction
//      )
//  {
//    Set<Set<Language>>  reconstruction = SymmetricDiff.clades(_reconstruction);
//    CounterMap<Set<Set<Language>>,PartialCoalescentState> organizedCounts
//      = new CounterMap<Set<Set<Language>>,PartialCoalescentState>();
//    for (PartialCoalescentState p : particles.keySet())
//      organizedCounts.incrementCount(p.allClades(), p, 
//          Purity.purity(Arbre.arbre2Tree(p.getUnlabeledArbre()), referenceClusters));
//    for (Set<Set<Language>> s : organizedCounts)
//      ;
//  }
  
//  private <S,T> Counter<S> marginals(CounterMap<S,T> cm)
//  {
//    Counter<S> result = new Counter<S>();
//    for (S key : cm.keySet())
//      result.setCount(key, cm.getCounter(key).totalCount());
//    return result;
//  }

//  private Coalescent _guideCoalescent = null;
  
//  public static PartialCoalescentState guidedPartialCoalescent(
//      String newickFile,
//      PartialCoalescentState initial)
//  {
//    return guidedPartialCoalescent(getGuideCoalescent(newickFile), initial);
//  }
  
//  public static PartialCoalescentState guidedPartialCoalescent(
//      final Coalescent guide, 
//      final PartialCoalescentState initial)
//  {
//    final Map<Set<Language>, Language> clades2internal = new HashMap<Set<Language>,Language>();
//    clades2internal(guide.topology(), clades2internal);
//    // find height of all the nodes
//    Arbre<Pair<Set<Language>,Double>> heights = guide.topology().postOrderMap(new ArbreMap<Language, Pair<Set<Language>,Double>>() {
//      @Override public Pair<Set<Language>, Double> map(Arbre<Language> currentDomainNode)
//      {
//        final Language curlang = currentDomainNode.getContents();
//        Set<Language> clade = new HashSet<Language>();
//        if (currentDomainNode.isLeaf())clade.add(curlang);
//        for (Pair<Set<Language>,Double> child : getChildImage())
//          clade.addAll(child.getFirst());
//        if (currentDomainNode.isLeaf()) return Pair.makePair(clade, 0.0);
//        Pair<Set<Language>,Double> firstChild = getChildImage().get(0);
//        return Pair.makePair(clade, firstChild.getSecond() 
//            + guide.branchLengths().get(clades2internal.get(firstChild.getFirst())));
//      }
//    });
//    // sotre them in a hash map
//    Map<Set<Language>,Double> clade2height = new HashMap<Set<Language>,Double>();
//    for (Pair<Set<Language>, Double> node : heights.nodeContents())
//      clade2height.put(node.getFirst(),node.getSecond());
//    // process all the coalescent events
//    PartialCoalescentState current = initial;
//    while (!current.isFinalState())
//    {
//      // figure out which pair to merge: the smallest
//      Counter<Pair<Integer,Integer>> candidates = new Counter<Pair<Integer,Integer>>();
//      for (int i = 0; i < current.nRoots(); i++)
//        for (int j = i+1; j < current.nRoots(); j++)
//        {
//          Set<Language> currentClade = new HashSet<Language>();
//          currentClade.addAll(current.clade(i));
//          currentClade.addAll(current.clade(j));
//          if (clade2height.containsKey(currentClade))
//            candidates.setCount(Pair.makePair(i,j), -clade2height.get(currentClade));
//        }
//      Pair<Integer,Integer> best = candidates.argMax();
//      // merge it
//      current = current.coalesce(best.getFirst(), best.getSecond(), -candidates.getCount(best));
//    }
//    return current;
//  }
//  
//  private static Set<Language> clades2internal(
//      Arbre<Language> topology, Map<Set<Language>,Language> map)
//  {
//    Set<Language> result = new HashSet<Language>();
//    if (topology.isLeaf()) result.add(topology.getContents());
//    else for (Arbre<Language> child : topology.getChildren())
//      result.addAll(clades2internal(child, map));
//    map.put(result, topology.getContents());
//    return result;
//  }

  public static RootedTree getGuideCoalescent(
      final String guideTreePath)
  {
    return getGuideCoalescent(guideTreePath, null,1.0);
  }
  public static RootedTree getGuideCoalescent(
      final String guideTreePath, 
      final Dataset data,
      final double guideTreeScaling)
  {
//    if (_guideCoalescent != null) return _guideCoalescent;
    RootedTree _guideCoalescent;
    if (guideTreePath.equals(autoGuideTree))
    {
      Map<Taxon,String> clusters = data.getReferenceClusters();
      final Map<Taxon,Double> bl = new HashMap<Taxon,Double>();
      List<Arbre<Taxon>> children = new ArrayList<Arbre<Taxon>>();
      Set<Taxon> obsLang = data.observations().keySet();
      for (String cluster : new HashSet<String>(clusters.values()))
      {
        List<Arbre<Taxon>> children2 = new ArrayList<Arbre<Taxon>>();
        for (Taxon lang : clusters.keySet())
          if (obsLang.contains(lang))
            if (clusters.get(lang).equals(cluster))
            {
              Arbre<Taxon> leaf = Arbre.arbre(lang);
              children2.add(leaf);
              bl.put(lang, 1.0);
            }
        if (children2.size() > 0)
        {
          Taxon curLang = new Taxon(cluster);
          children.add(Arbre.arbre(curLang, children2));
          bl.put(curLang, 2.0);
        }
      }
      final Arbre<Taxon> root = Arbre.arbre(new Taxon("root"), children);
      _guideCoalescent = new RootedTree() {
        public Arbre<Taxon> topology()            { return root;  }
        public Map<Taxon, Double> branchLengths() { return bl; } 
        @Override
        public int nTaxa()
        {
          return root.nLeaves();
        }
        @Override
        public RootedTree getRooted()
        {
          return this;
        }
        @Override
        public UnrootedTree getUnrooted()
        {
          return UnrootedTree.fromRooted(this);
        }
      
      } ;
    }
    else try
    {
      final RootedTree c = RootedTree.Util.load(new File(guideTreePath));
//      NewickParser np = new NewickParser(IOUtils.openIn(guideTreePath));
//      Tree<String> tree = np.parse();
//      final Arbre<Language> a = Arbre.tree2Arbre(tree).preOrderMap(
//          new ArbreMap<String,Language>() {
//            @Override  public Language map(Arbre<String> c) 
//            { return new Language(c.getContents()); } });
      final Map<Taxon,Double> bl = new HashMap<Taxon,Double>();
      for (Taxon lang : c.branchLengths().keySet())
        bl.put(lang, c.branchLengths().get(lang) * guideTreeScaling);
      _guideCoalescent = new RootedTree() {
        public Arbre<Taxon> topology()            { return c.topology();  }
        public Map<Taxon, Double> branchLengths() { return bl; }
        @Override
        public int nTaxa()
        {
          return c.topology().nLeaves();
        }
        @Override
        public RootedTree getRooted()
        {
          return this;
        }
        @Override
        public UnrootedTree getUnrooted()
        {
          return UnrootedTree.fromRooted(this);
        }
      };
    } catch (Exception e) { throw new RuntimeException(e); }
    outputTree(_guideCoalescent.topology(), "guide-tree");
    return _guideCoalescent;
  }
  public static Arbre<Taxon> outputTree(Arbre<Taxon> reconstruction, String prefix)
  {
    return outputTree(reconstruction, prefix, null);
  }
  public static Arbre<Taxon> outputTree(Arbre<Taxon> reconstruction, String prefix, Map<Taxon,Double> bl)
  {
    IO.writeToDisk(
        Execution.getFile(prefix + ".newick"), 
        DataPrepUtils.newick(Arbre.arbre2Tree(reconstruction),bl,true));
    IO.writeToDisk(
        Execution.getFile(prefix + ".txt"), 
        reconstruction.deepToString());
    LogInfo.logs("Reconstructed tree (" + prefix + "):\n" + reconstruction.deepToString());
    return reconstruction;
  }
  
  public static void main(String [] args)
  { 
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    if (!Arrays.asList(args).contains("NOJARS"))
      Execution.jarFiles = new ArrayList<String>(Arrays.asList(
          "/home/eecs/bouchard/jars/ptychodus.jar",
          "/home/eecs/bouchard/jars/nuts.jar",
          "/home/eecs/bouchard/jars/pepper.jar",
          "/home/eecs/bouchard/jars/fig.jar"));
    pf = new ParticleFilter<PartialCoalescentState>();
    loader = new CTMCLoader();
    Execution.run(args, new Train(), 
        "wals", WalsDataset.class, 
        "hgdb", HGDPDataset.class,
        "dataGen", GeneratedDataset.class,
        "paramGen", GeneratedDataset.genCTMCLoader,
        "filter", pf, 
        "init", loader,
        "suffstat", LearningProcessor.class);
  }
}

//public class Train implements Runnable
//{
//  @Option public static int nEMIters = 20;
//  @Option public static DatasetType datasetType = DatasetType.WALS;
//  @Option public static boolean useReversibleMtx = true;
//  @Option public static int increaseNSamplesPerEMIteration = 0;
//  @Option public static EstimationMethod estimationMethod = EstimationMethod.FULL;
//  
//  public static enum EstimationMethod 
//  {
//    DONT {
//      @Override
//      public CTMC estimate(CTMC old, LearningProcessor learningProcessor, int emIter) { return old; }
//      @Override
//      public boolean learningProcessorNeeded()
//      {
//        // TODO Auto-generated method stub
//        return false;
//      }
//    }, 
//    FULL {
//      @Override
//      public CTMC estimate(CTMC old, LearningProcessor learningProcessor, int emIter)
//      {
//        LogInfo.track("Reestimating parameters");
//        double [] sd = old.getInitialDistribution(0);
//        if (Train.useReversibleMtx) LogInfo.logs("Stationary distribution:" + Arrays.toString(sd));
//        CTMC result = (Train.useReversibleMtx ? learningProcessor.getTiedReversibleMLE(sd) : learningProcessor.getTiedMLE());
//        CTMCUtils.saveInExec(result, "reest-" + (emIter+1));
//        LogInfo.end_track();
//        return result;
//      }
//    }, CONSTRAINED;
//    abstract CTMC estimate(CTMC old, LearningProcessor learningProcessor, int emIter);
//    abstract boolean learningProcessorNeeded();
//  }
//  
//  private static ParticleFilter<PartialCoalescentState> pf;
//  private static CTMCLoader loader;
//  
//  private CTMC ctmc; // current params
//  
//  public void run()
//  {
//    // load data, init params
//    Dataset data = datasetType.loadDataset(); 
//    loader.setData(data);
//    ctmc = loader.load();
//    CTMCUtils.saveInExec(ctmc, "init");
//    // do EM
//    for (int emIter = 0; emIter < nEMIters; emIter++)
//    {
//      LogInfo.track("EM Iteration " + (emIter+1) + "/" + nEMIters,true);
//      // prepare sampling machinery
//      PartialCoalescentState initState = PartialCoalescentState.initState(data, ctmc);
//      ParticleMapperProcessor<PartialCoalescentState, Set<Set<Language>>> mbr =
//        SymmetricDiff.createCladeProcessor();
//      LearningProcessor ssp = null;
//      ForkedProcessor<PartialCoalescentState> processors 
//        = new ForkedProcessor<PartialCoalescentState>(mbr);
//      if (estimationMethod.learningProcessorNeeded())
//      {
//        ssp = new LearningProcessor(ctmc);
//        processors.processors.add(ssp);
//      }
//      PriorPriorKernel kernel = new PriorPriorKernel(initState);
//      // E step : sample!
//      pf.bootstrapFilter(kernel, processors);
//      pf.N += increaseNSamplesPerEMIteration;
//      // reconstruct a tree
//      Arbre<Language> reconstruction = outputTree(mbr.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE), "consensusTree-" + emIter);
//      // evaluate purity if possible
//      if (data.hasReferenceClusters())
//        LogInfo.logs("Purity:" + 
//            Purity.purity(Arbre.arbre2Tree(reconstruction), 
//                          data.getReferenceClusters()));
//      // M step : update parameters
//      ctmc = estimationMethod.estimate(this);
//      LogInfo.end_track();
//    }
//  }
//  
//  public static Arbre<Language> outputTree(Set<Set<Language>> clades, String prefix)
//  {
//    Arbre<Language> reconstruction = SymmetricDiff.clades2arbre(clades);
//    IO.writeToDisk(
//        Execution.getFile(prefix + ".newick"), 
//        DataPrepUtils.newick(Arbre.arbre2Tree(reconstruction)));
//    IO.writeToDisk(
//        Execution.getFile(prefix + ".txt"), 
//        reconstruction.deepToString());
//    LogInfo.logs("Reconstructed tree (" + prefix + "):\n" + reconstruction.deepToString());
//    return reconstruction;
//  }
//  
//  public static void main(String [] args)
//  { 
//    Execution.monitor = true;
//    Execution.makeThunk = false;
//    Execution.create = true;
//    Execution.useStandardExecPoolDirStrategy = true;
//    if (!Arrays.asList(args).contains("NOJARS"))
//      Execution.jarFiles = new ArrayList<String>(Arrays.asList(
//          "/home/eecs/bouchard/jars/ptychodus.jar",
//          "/home/eecs/bouchard/jars/nuts.jar",
//          "/home/eecs/bouchard/jars/pepper.jar",
//          "/home/eecs/bouchard/jars/fig.jar"));
//    pf = new ParticleFilter<PartialCoalescentState>();
//    loader = new CTMCLoader();
//    Execution.run(args, new Train(),
//        "train", Train.class,
//        "wals", WalsDataset.class, 
//        "hgdb", HGDPDataset.class,
//        "dataGen", GeneratedDataset.class,
//        "paramGen", GeneratedDataset.genCTMCLoader,
//        "filter", pf, 
//        "init", loader,
//        "suffstat", LearningProcessor.class);
//  }
//}

