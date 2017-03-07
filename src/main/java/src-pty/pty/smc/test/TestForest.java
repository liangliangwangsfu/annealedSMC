package pty.smc.test;
import java.io.*;
import java.util.*;

import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Hasher;

import pty.Train;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.WalsData;
import pty.io.WalsDataset;
import pty.io.WalsDataset.LanguageDatabase;
import pty.learn.CTMCLoader;
import pty.smc.ConstrainedKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.ForestModelCalculator;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.Taxon;

public class TestForest  implements Runnable
{
  @Option public double langInvRate = 10;
  @Option public double rootHeight = 10.0;
  private static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
  private static CTMCLoader langParamLoader = new CTMCLoader();
  private Dataset langData;
  
  public static void main (String args []){
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    Execution.run(args, new TestForest(),
        "wals", WalsDataset.class,
        "kernel", ConstrainedKernel.class,
        "filter", pf,
        "pcs", PartialCoalescentState.class,
        "ppk", PriorPriorKernel.class);
  }
  
  
  private PartialCoalescentState initLanguageState () {
    langData = WalsDataset.getPreprocessedCorpus();
    langParamLoader.setData(langData);
    CTMC langParam = langParamLoader.load();
    return PartialCoalescentState.initForestState(langData, langParam, rootHeight, langInvRate); 
  }

  @SuppressWarnings("unchecked")
  public void run()
  {
    PartialCoalescentState initState = initLanguageState();
    ParticleKernel<PartialCoalescentState> ppk = new ConstrainedKernel(initState, getConstraints());
    SingleTreePosteriorDecoder dec = new SingleTreePosteriorDecoder();
    ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> mbr =
      SymmetricDiff.createCladeProcessor();
    ForkedProcessor<PartialCoalescentState> processors 
    = new ForkedProcessor<PartialCoalescentState>(mbr,dec);
    pf.sample ( ppk , processors);
    Train.outputTree( 
        SymmetricDiff.clades2arbre(
            mbr.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE)), 
            "consensusTree");
    LogInfo.logs("Single tree posterior:" + dec.getSingleTreePosterior());
  }
  
  private Set<Set<Taxon>> getConstraints()
  {
    Set<Set<Taxon>> constraints = new HashSet<Set<Taxon>>();
    constraints.addAll(CollUtils.inducedPartition(WalsDataset.langDB.genusMap()));
    constraints.addAll(CollUtils.inducedPartition(WalsDataset.langDB.familyMap()));
    LogInfo.logs("Constraints:" + constraints);
    return constraints;
  }
  
  public static class SingleTreePosteriorDecoder implements ParticleProcessor<PartialCoalescentState>
  {
    private double singleTreePosterior = 0.0; 
    public void process(PartialCoalescentState state, double weight)
    {
      ForestModelCalculator node = (ForestModelCalculator) (state.getLikelihoodModelCalculator(0));
      singleTreePosterior += weight * (1.0 - node.posteriorNoLanguagePr());
    }
    public double getSingleTreePosterior() { return singleTreePosterior; }
  }
}
