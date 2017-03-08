package pty;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.omg.CORBA.MARSHAL;

import nuts.io.IO;
import nuts.lang.ArrayUtils;
import nuts.math.GMFct;
import nuts.util.Arbre;
import pty.eval.Purity;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.learn.CTMCLoader;
import pty.learn.DiscreteBP;
import pty.learn.LearningProcessor;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.MAPDecoder;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.CTMCUtils;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.StrUtils;
import fig.exec.Execution;
import goblin.DataPrepUtils;
import goblin.Taxon;

public class ProcessPops implements Runnable
{
  private static ParticleFilter<PartialCoalescentState> pf;
  private static CTMCLoader loader;
  
  
  @SuppressWarnings("unchecked")
  public void run()
  {
    Map<Taxon,double[]> result = new HashMap<Taxon,double[]>();
    // load data, init params
    final Dataset fullData = DatasetType.HGDP.loadDataset(); 
    loader.setData(fullData);
    CTMC ctmc = loader.load();
    CTMCUtils.saveInExec(ctmc, "init");
    final Map<Taxon,String> pops = fullData.getReferenceClusters();
    Set<String> allPopulations = new HashSet<String>(pops.values());
    Map<Taxon,double[][]> allData = fullData.observations();
    int i = 1;
    for (String currentPopulation : allPopulations)
    {
      LogInfo.track("Processing population " + currentPopulation + " (" + (i++) + "/" + allPopulations.size() + ")",true);
      // construct a mini-dataset containing only the individual in the current population
      final Map<Taxon,double[][]> _restrictedData = dataForPopulation(pops, allData, currentPopulation);
      Dataset restrictedData = new Dataset() {
        public Map<Taxon, String> getReferenceClusters() { return pops; }
        public boolean hasReferenceClusters() { return true; }
        public Map<Taxon, double[][]> observations() { return _restrictedData; }
        public int nCharacter(int site) { return fullData.nCharacter(site); }
        public int nSites() { return fullData.nSites(); }
      };
      PartialCoalescentState initState = PartialCoalescentState.initState(restrictedData, ctmc);
      // prepare sampling machinery
      ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> mbr =
        SymmetricDiff.createCladeProcessor();
      MAPDecoder<PartialCoalescentState> mapDecoder = new MAPDecoder<PartialCoalescentState>(); 
      ForkedProcessor<PartialCoalescentState> processors = new ForkedProcessor<PartialCoalescentState>(mbr,mapDecoder);
      PriorPriorKernel kernel = new PriorPriorKernel(initState);
      // E step : sample!
      pf.sample(kernel, processors);
      // add a new pop in the processed pop file
      result.put(new Taxon(currentPopulation), getRootPosterior(mapDecoder.map()));
      // reconstruct a tree for fun
      Train.outputTree( SymmetricDiff.clades2arbre(mbr.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE)), "consensusTree-" + currentPopulation);
      LogInfo.end_track();
    }
    writeResult(result);
  }
  

  /**
   * Output is similar to PHYLIP format, except that real entries are used to 
   * represent the posterior on 2 sites
   * @param result
   */
  private void writeResult(Map<Taxon, double[]> result)
  {
    String file = Execution.getFile("result");
    PrintWriter out = IOUtils.openOutHard(file);
    for (Taxon lang : result.keySet())
    {
      out.append(lang + "\t");
      out.append(StrUtils.join(result.get(lang), "\t") + "\n");
    }
    out.close();
  }



  private double [] getRootPosterior(PartialCoalescentState state)
  {
    final int nSites = state.getCTMC().nSites();
    Arbre<Taxon> a = state.getArbreAndBranchLengths().getFirst();
    Taxon root = a.getContents();
    double[] result = new double[nSites];
    for (int s = 0; s < nSites; s++)
    {
      GMFct<Taxon> posterior = DiscreteBP.posteriorMarginalTransitions(state, s);
      result[s] = posterior.get(root, 0);
    }
    return result;
  }



  private Map<Taxon, double[][]> dataForPopulation(Map<Taxon, String> pops,
      Map<Taxon, double[][]> allData, String currentPopulation)
  {
    Map<Taxon,double[][]> result = new HashMap<Taxon,double[][]>();
    for (Taxon lang : allData.keySet())
      if (pops.get(lang).equals(currentPopulation))
        result.put(lang,allData.get(lang));
    return result;
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
    Execution.run(args, new ProcessPops(), 
        "hgdb", HGDPDataset.class,
        "filter", pf, 
        "init", loader);
  }
}
