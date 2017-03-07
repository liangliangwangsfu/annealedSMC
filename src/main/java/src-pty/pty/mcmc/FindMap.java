package pty.mcmc;
import java.io.*;
import java.util.*;

import pty.UnrootedTree;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.learn.CTMCLoader;
import pty.mcmc.PhyloSampler.NonClockTreePrior;
import pty.smc.MapLeaves;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.models.BrownianModel;
import pty.smc.models.CTMC;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Parallelizer;
import fig.exec.Execution;
import goblin.Taxon;
import goblin.BayesRiskMinimizer.LossFct;

import nuts.io.IO;
import nuts.util.CollUtils;

public class FindMap implements Runnable
{
  @Option(required=true) public String experimentDir = null;
  @Option public int numberOfThreads = 1;
  
  private static Main main = new Main();
  private MapLeaves ml = null;
  private static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
  private static CTMCLoader loader = new CTMCLoader();  
  private Dataset wals, hgdp;
  
  /**
   * @param args
   */
  public static void main(String[] args)
  {
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    Execution.run(args, new FindMap(),
        "main", main,
        "prior", PhyloSampler._defaultPriorOptions,
        "prop", ProposalDistribution.Util._defaultProposalDistributionOptions,
        "sampler", PhyloSampler._defaultPhyloSamplerOptions,
        "wals", WalsDataset.class, 
        "hddp", HGDPDataset.class,
        "filter", pf,
        "langparam", loader,
        "partemp", ParallelTemperingChain._defaultTemperingOptions);
  }

  private static Set<String> fileNames(List<File> subFiles)
  {
    Set<String> result = CollUtils.set();
    for (File f : subFiles)
      result.add(f.getName());
    return result;
  }
  
  private double max = Double.NEGATIVE_INFINITY;
  private UnrootedTree
    bestWalsT = null,
    bestHgdpT = null;

  @Override
  public void run()
  {
    // load trees
    List<File> 
      subFiles1 =  IO.locate(new File(experimentDir, "WALS-chain-0"), IO.suffixFilter("gz","newick")),
      subFiles2 =  IO.locate(new File(experimentDir, "HGDP-chain-0"), IO.suffixFilter("gz","newick"));
    final Set<String> fileNames = fileNames(subFiles1);
    fileNames.retainAll(fileNames(subFiles2));
    // data
    wals = DatasetType.WALS.loadDataset();
    hgdp = DatasetType.HGDP.loadDataset();
    loader.setData(wals);
    final CTMC ctmc = loader.load();
    final BrownianModel bm = new BrownianModel(hgdp.nSites(), main.brownianMotionVariance);
    ml = MapLeaves.parse (main.mapfile);
    // find map
    final NonClockTreePrior prior = PhyloSampler._defaultPhyloSamplerOptions.prior.prior(PhyloSampler._defaultPriorOptions);
    LogInfo.track("Searching map");
    Parallelizer<String> parallelizer = new Parallelizer<String>(numberOfThreads);
    parallelizer.setPrimaryThread();
    parallelizer.process(new ArrayList<String>(fileNames), new Parallelizer.Processor<String>() {
      public void process(String name, int _i, int _n, boolean log) {
        LogInfo.logsForce("Processed file " + (_i+1) + "/" + fileNames.size());
        File 
          walsFile = new File(experimentDir + "/WALS-chain-0/" + name),
          hgdpFile = new File(experimentDir + "/HGDP-chain-0/" + name);
        try
        {
          for (List<String> lines : IO.i(walsFile, hgdpFile))
          {
            // create trees
            String 
              walsLine = lines.get(0),
              hgdpLine = lines.get(1);
            UnrootedTree 
              walsTree = UnrootedTree.fromNewick(walsLine),
              hgdpTree = UnrootedTree.fromNewick(hgdpLine);
            // create likelihood models
            UnrootedTreeState 
              walsState = UnrootedTreeState.fromCTMC(walsTree, wals, ctmc),
              hgdpState = UnrootedTreeState.fromBrownianMotion(hgdpTree, hgdp, bm);
            double sum = walsState.logLikelihood() + hgdpState.logLikelihood();
            // prior (base measure)
            sum += prior.logPriorDensity(walsTree) + prior.logPriorDensity(hgdpTree);
            // prior (agreement)
            sum += - main.initialAgreementParam 
              * SymmetricDiff.symmetricDifferenceSize(
                  ml.filterClades(walsTree.clades()), 
                  ml.mapClades(ml.filterClades(hgdpTree.clades())));
            synchronized (main)
            {
              if (sum > max)
              {
                LogInfo.logs("Higher score found: " + sum);
                max = sum;
                bestWalsT = walsTree;
                bestHgdpT = hgdpTree;
              }
            }
          }
        }
        catch (Exception e) { LogInfo.warning("Problem with file " + name + "... skipping"); }
      }});
    LogInfo.logsForce("Best score: " + max);
    LogInfo.logsForce("Best wals: " + bestWalsT.toNewick());
    LogInfo.logsForce("Best hgdp: " + bestHgdpT.toNewick());
  }

}
