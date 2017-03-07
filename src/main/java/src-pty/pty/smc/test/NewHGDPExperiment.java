package pty.smc.test;
import java.io.*;
import java.util.*;

import ev.ex.TreeGenerators;
import ev.poi.processors.TreeDistancesProcessor;
import ev.to.NJ;
import fig.basic.Option;
import fig.exec.Execution;
import gep.util.OutputManager;
import goblin.Taxon;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.Dataset.DatasetType;
import pty.mcmc.UnrootedTreeState;
import pty.mcmc.ParallelTemperingChain;
import pty.mcmc.PhyloSampler;
import pty.mcmc.PhyloSampler.PhyloProcessor;
import pty.smc.ConstrainedKernel;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.models.BrownianModel;
import pty.smc.models.BrownianModelCalculator;

import nuts.io.IO;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class NewHGDPExperiment implements Runnable
{
  @Option public HGDPExperimentType type = HGDPExperimentType.SMC_NONCLOCK;
  @Option public Random inferenceRand = new Random(1);
  @Option public int nIterations = 1000;
  @Option public int parallelism = 1;
  @Option public int trivialParallelism = 1;
  @Option public boolean resampleRoot = false;
  
  public static enum HGDPExperimentType
  {
    MCMC_NONCLOCK {

      @Override
      public void doIt(NewHGDPExperiment xp)
      {
        // random init tree
        PartialCoalescentState state = PartialCoalescentState.initState(xp.data, xp.bm, false);
        List<Taxon> leaves = list(state.getObservations().observations().keySet());
        RootedTree initTree = TreeGenerators.sampleExpNonclock(xp.inferenceRand, leaves, NCPriorPriorKernel.deltaProposalRate);
        UnrootedTreeState ncts = UnrootedTreeState.fromBrownianMotion(UnrootedTree.fromRooted(initTree), xp.data, xp.bm);
        
        PhyloSampler sampler = new PhyloSampler();
        sampler.init(ncts);
        
        PhyloProcessor tdpa = 
          new PhyloSampler.PhyloProcessorAdaptor(xp.processor);
        sampler.setProcessors(Arrays.asList(tdpa));
        
        ParallelTemperingChain temperingChain = new ParallelTemperingChain();
        temperingChain.init(sampler);
        temperingChain.sample();
      }
      
    },
    SMC_NONCLOCK {
      @Override
      public void doIt(NewHGDPExperiment xp)
      {
        PartialCoalescentState state = PartialCoalescentState.initState(xp.data, xp.bm, xp.resampleRoot);
        ParticleKernel<PartialCoalescentState> ppk = new NCPriorPriorKernel(state);
        TestBrownianModel.pf.sample(ppk,xp.processor);
      }
    };
    public abstract void doIt(NewHGDPExperiment xp);
  }
  
  public static void main(String [] args)
  {
    IO.run(args, new NewHGDPExperiment(),        
        "hgdp", HGDPDataset.class,
        "kernel", ConstrainedKernel.class,
        "filter", TestBrownianModel.pf,
        "bm", BrownianModelCalculator.class,
        "pcs", PartialCoalescentState.class,
        "nc", NCPriorPriorKernel.class,
        "nj", NJ.class);
  }
  
  private OutputManager om = new OutputManager();
//  private PartialCoalescentState geneState;
//  private ParticleKernel<PartialCoalescentState> ppk;
  private TreeDistancesProcessor processor;
  private Dataset data;
  private BrownianModel bm;

  @Override
  public void run()
  {
    TestBrownianModel.pf.rand = inferenceRand;
    TestBrownianModel.pf.nThreads = parallelism;
    TestBrownianModel.pf.N = nIterations;
//    PhyloSampler._defaultPhyloSamplerOptions.nIteration = nIterations;
    PhyloSampler._defaultPhyloSamplerOptions.rand = inferenceRand;
    ParallelTemperingChain._defaultTemperingOptions.rand = inferenceRand;
    ParallelTemperingChain._defaultTemperingOptions.nSwapsPerRoundPerChain = 1;
    ParallelTemperingChain._defaultTemperingOptions.nChainMovesPerRoundPerChain = 10;
    ParallelTemperingChain._defaultTemperingOptions.nRounds = nIterations / parallelism / ParallelTemperingChain._defaultTemperingOptions.nChainMovesPerRoundPerChain;
    ParallelTemperingChain._defaultTemperingOptions.nThreads = parallelism;
    ParallelTemperingChain._defaultTemperingOptions.nChains = parallelism;
    data = DatasetType.HGDP.loadDataset();
    bm = new BrownianModel(data.nSites(), 1.0);
    
    processor = new TreeDistancesProcessor();
    
    // do not actually parallelize this for now (would need to create copies of all the objects, etc)
    final long beg = System.currentTimeMillis();
    for (int i = 0; i < trivialParallelism; i++)
      type.doIt(this);
    final long delta = System.currentTimeMillis() - beg;
    
    UnrootedTree inferred = processor.getConsensus();
    outputTree(inferred, "", delta);
    
//    UnrootedTree mapTree = processor.getMode().getPhylogeny().getUnrooted();
    UnrootedTree mapTree = processor.getMode().getUnrooted();
      //UnrootedTree.fromRooted(processor.get.getFullCoalescentState());
    outputTree(mapTree, "-map", delta);
    
    om.close();
  }
  
  private void outputTree(UnrootedTree inferred, String suffix, long delta)
  {
    IO.writeToDisk(new File(Execution.getFile("inferred" + suffix + ".newick")), inferred.toNewick());
    // TODO evaluate the likelihood of the inferred tree
    UnrootedTreeState ncs = UnrootedTreeState.fromBrownianMotion(inferred, data, bm);
    om.printWrite("likelihood" + suffix , "time", delta, "likelihood" + suffix, ncs.logLikelihood());
  }

}
