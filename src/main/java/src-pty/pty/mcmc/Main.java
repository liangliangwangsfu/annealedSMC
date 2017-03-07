package pty.mcmc;
import java.io.File;
import java.sql.DataTruncation;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import nuts.math.Sampling;
import nuts.util.CollUtils;
import nuts.util.EasyFormat;

import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import pty.UnrootedTree;
import pty.eval.CladeOverlap;
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
import pty.smc.ParticleKernel;
import pty.smc.models.BrownianModel;
import pty.smc.test.TestBrownianModel.KernelType;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import fig.prob.Gamma;
import goblin.Taxon;
import goblin.BayesRiskMinimizer.LossFct;

public class Main implements Runnable
{
  @Option public ArrayList<DatasetType> dataTypes = new ArrayList<DatasetType>(Arrays.asList(DatasetType.WALS));
  @Option public KernelType kernelType = KernelType.PRIOR_PRIOR;
  @Option public double brownianMotionVariance = 0.2;
  @Option public ArrayList<String> pathsToInitTree = new ArrayList<String>(Arrays.asList(""));
  @Option public int nStepsPerSubRound = 150;
  @Option public Random paramRand = new Random(1);
  @Option public int nParamResamplingPerRound = 10;
  @Option public String mapfile = "data/language-gene-map.txt";
  @Option public boolean resampleAgreementParam = true;
  @Option public double paramResamplingStepSize = 2.0;
  @Option public double hyperParamScale = 3.0;
  @Option public double hyperParamShape = 2.0;
  @Option public double initialAgreementParam = 1;
  @Option public boolean monotonicallyIncreaseAgreementParam = false;
  @Option public double monotonicIncrease = 1;
  @Option public AgreementLoss cladeLoss = AgreementLoss.CLADEOVER;
  
  @Option public boolean conditionOnFamilies = true;
//  @Option public boolean conditionOnGenera = true;
  
  private MapLeaves ml = null;
  private static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
  private static CTMCLoader loader = new CTMCLoader();

  private Dataset data;
  
  public static void main (String args [])
  {
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
    Execution.run(args, new Main(),
        "prior", PhyloSampler._defaultPriorOptions,
        "prop", ProposalDistribution.Util._defaultProposalDistributionOptions,
        "sampler", PhyloSampler._defaultPhyloSamplerOptions,
        "wals", WalsDataset.class, 
        "hddp", HGDPDataset.class,
        "filter", pf,
        "langparam", loader,
        "partemp", ParallelTemperingChain._defaultTemperingOptions);
  }
  private Set<Taxon> allLangs = new HashSet<Taxon>();
  public void run() { run(""); }
  /**
   * return File objects pointing to the output samples
   * @param prefixForSampleDirectories
   * @return
   */
  public List<File> run(String prefixForSampleDirectories)
  {
    List<File> result = new ArrayList<File>(); 
    Map<DatasetType, ParallelTemperingChain> chains = CollUtils.map();
    for (int i = 0; i < dataTypes.size(); i++)
    {
      DatasetType dataType = dataTypes.get(i);
      String pathToInitTree = "";
      if (pathsToInitTree.size() > i)
        pathToInitTree = pathsToInitTree.get(i);
      data = dataType.loadDataset(); 
      loader.setData(data);
      
      allLangs.addAll(data.observations().keySet());
      
      UnrootedTreeState nctsInit = null;
      
      if (pathToInitTree.equals(""))
      {
        // run the pf
        LogInfo.track("Particle Filter phase");
        PartialCoalescentState initState = null;
        if (dataType == DatasetType.HGDP)
          initState = PartialCoalescentState.initState(data, new BrownianModel(data.nSites(), brownianMotionVariance), false);
        else if (dataType == DatasetType.WALS)
          initState = PartialCoalescentState.initState(data, loader.load());
        else throw new RuntimeException();
        PartialCoalescentState initTree = null;
        ParticleKernel<PartialCoalescentState> pk = kernelType.load(initState, null);
        ParticleFilter.MAPDecoder<PartialCoalescentState> processor 
          = new ParticleFilter.MAPDecoder<PartialCoalescentState>();
        pf.sample( pk, processor);
        initTree = processor.map();
        LogInfo.end_track();
        nctsInit = UnrootedTreeState.fromPartialCoalescentState(initTree);
      }
      else
      {
        UnrootedTree nct = UnrootedTree.fromNewick(new File(pathToInitTree));
        if (dataType == DatasetType.HGDP)
          nctsInit = UnrootedTreeState.fromBrownianMotion(nct, data, new BrownianModel(data.nSites(), brownianMotionVariance));
        else if (dataType == DatasetType.WALS)
          nctsInit = UnrootedTreeState.fromCTMC(nct, data, loader.load());
        else throw new RuntimeException();
      }
      
      // start MCMC
      
      PhyloSampler sampler = new PhyloSampler();
      sampler.init(nctsInit);
      
      if (dataType == DatasetType.WALS)
      {
        if (conditionOnFamilies)
          sampler.setConstraints(WalsDataset.langDB.familyMap());
//        if (conditionOnGenera)
//          sampler.addConstraints(WalsDataset.langDB.genusMap());
      }
      
      ParallelTemperingChain temperingChain = new ParallelTemperingChain();
      result.add(temperingChain.setOutputPrefix(prefixForSampleDirectories + "-" + dataType + "-"));
      temperingChain.init(sampler);
      chains.put(dataType, temperingChain);
    }
    if (chains.size() == 1)
    {
      LogInfo.track("Markov Chain Monte Carlo phase");
      chains.values().iterator().next().sample();
      LogInfo.end_track();
    }
    else if (chains.size() > 2) throw new RuntimeException();
    else
      sampleJoint(chains);
    return result;
  }

  private int nOuterRuns = -1;

  private void sampleJoint(Map<DatasetType, ParallelTemperingChain> chains)
  {
    ml = MapLeaves.parse (mapfile);
    if (monotonicallyIncreaseAgreementParam && resampleAgreementParam)
      throw new RuntimeException();
    
//    {
//      System.out.println("HACK!");
//      DatasetType src = DatasetType.HGDP, dst = DatasetType.WALS;
//      // take the src tree
//      NonClockTree nct = chains.get(src).getChain(0).getInitialState().getNonClockTree();
//      // relabel the leaves
//      
//    }
    
    LogInfo.track("Checking map file",true);
    LogInfo.logs("Data not in map items: " + ml.dataNotInMapItems(allLangs));
    LogInfo.logs("Map items not in data: " + ml.mapItemsNotInData(allLangs));
    LogInfo.end_track();
    
    List<JointPrior<Set<Set<Taxon>>>> jointPriors = initJoints(chains);
    // iterate b/w sampling, updating joint prior, and resampling agreement parameters
    update(jointPriors, chains); // finished initialization!
    for (int iter = 0; iter < nOuterRuns; iter++)
      for (DatasetType dt : chains.keySet())
      {
        chains.get(dt).sample();
        update(jointPriors, chains);
//        LogInfo.logsForce("Current agreement: " + jointPriors.get(0).getCurrentDistance());
        LogInfo.track("Agreement statistics:");
        LogInfo.logsForce(toString(jointPriors));
        LogInfo.end_track();
        if (resampleAgreementParam || monotonicallyIncreaseAgreementParam)
          resampleParams(jointPriors);
      }
  }
  
  
  
  public String toString(List<JointPrior<Set<Set<Taxon>>>> jointPriors)
  {
    StringBuilder result = new StringBuilder();
    for (int i = 0; i < jointPriors.size(); i++)
      result.append("#"  + i + ": " +
          "agreeParam=" + EasyFormat.fmt2(jointPriors.get(i).agreementParameter) +", " + 
          "dist=" + EasyFormat.fmt2(jointPriors.get(i).getCurrentDistance()) + ", " +
          "paramAccept=" + EasyFormat.fmt2(jointPriors.get(i).resamplingStats.getMean()) 
          + "\n");
    return result.toString();
  }
  
  private void resampleParams(List<JointPrior<Set<Set<Taxon>>>> jointPriors)
  {
    for (JointPrior<Set<Set<Taxon>>> jp : jointPriors)
      if (monotonicallyIncreaseAgreementParam)
        jp.increaseAgreementParam(monotonicIncrease);
      else
        for (int iter = 0; iter < nParamResamplingPerRound; iter++)
          jp.sampleAgreementParameter(paramRand);
//    LogInfo.logsForce("Resampled param, accept rate: " + paramResamplingStats.getMean() + ", current value in main chain: " + 
//        jointPriors.get(0).agreementParameter);
  }

  private void update(List<JointPrior<Set<Set<Taxon>>>> jointPriors,
      Map<DatasetType, ParallelTemperingChain> chains)
  {
    for (int c = 0; c < chains.values().iterator().next().nChains(); c++)
    {
      JointPrior<Set<Set<Taxon>>> jp = jointPriors.get(c);
      for (DatasetType dt : chains.keySet())
        jp.update(dt, chains.get(dt).getChain(c).getCurrentState().getNonClockTree());
    }
  }
  
  public enum AgreementLoss 
  { 
    SYMMDIFF {
      @Override
      public LossFct<Set<Set<Taxon>>> getLoss(Main main)
      {
        return main.new SymmDiffLoss();
      }
    }, CLADEOVER {
      @Override
      public LossFct<Set<Set<Taxon>>> getLoss(Main main)
      {
        return main.new CladeOverlapLoss();
      }
    };
    abstract LossFct<Set<Set<Taxon>>> getLoss(Main main);
  }
  
  class SymmDiffLoss implements LossFct<Set<Set<Taxon>>>
  {
    @Override
    public double loss(Set<Set<Taxon>> t1, Set<Set<Taxon>> t2)
    {
      final double result =  SymmetricDiff.symmetricDifferenceSize(t1, ml.mapClades(t2)) / 2.0 / ((double)t1.size());
      return result;
    }
  }
  
  class CladeOverlapLoss implements LossFct<Set<Set<Taxon>>>
  {
    @Override
    public double loss(Set<Set<Taxon>> t1, Set<Set<Taxon>> t2)
    {
      final double n = t1.size();
      return 1.0 - CladeOverlap.cladeOverlap(t1, ml.mapClades(t2)) / n /n ;
    }
  }

  private List<JointPrior<Set<Set<Taxon>>>> initJoints(Map<DatasetType, ParallelTemperingChain> chains)
  {
    List<JointPrior<Set<Set<Taxon>>>> jointPriors = CollUtils.list(); // one for each chain
    for (int i = 0; i < chains.values().iterator().next().nChains(); i++)
    {
      JointPrior<Set<Set<Taxon>>> jp = new JointPrior<Set<Set<Taxon>>>(
          new CladeExtractor(),
          cladeLoss.getLoss(this),
          paramResamplingStepSize,
          hyperParamShape, hyperParamScale,
          initialAgreementParam);
      for (DatasetType dt : chains.keySet())
        jp.addToJointModel(dt, chains.get(dt).getChain(i));
      jointPriors.add(jp);
    }
    for (DatasetType dt : chains.keySet())
    {
      int initialNumberOfRounds = chains.get(dt).getOptions().nRounds; // use nStepsPerSubRound to change chain's number of round
      if (nOuterRuns == -1)
        nOuterRuns = initialNumberOfRounds/nStepsPerSubRound;
      chains.get(dt).getOptions().nRounds = nStepsPerSubRound;
    }
    return jointPriors;
  }
  
  
  public static class JointPrior<SuffStat>
  {
    public final NonClockTreeSuffStatExtractor<SuffStat> suffStatExtractor;
    public final LossFct<SuffStat> lossFct;
    public final double a;
    public final double hyperShape,hyperScale;
    private double temp = Double.NaN;
    private double agreementParameter;
    private Map<DatasetType, SuffStat> statesSS = CollUtils.map();
    private Map<DatasetType, UnrootedTree> states = CollUtils.map();
    private Map<DatasetType, NonClockTreePrior> basePriors = CollUtils.map();
    private SummaryStatistics resamplingStats = new SummaryStatistics();
//    private Map<DatasetType, PhyloSampler> baseSamplers = CollUtils.emMap();
    public JointPrior(
        NonClockTreeSuffStatExtractor<SuffStat> suffStatExtractor,
        LossFct<SuffStat> lossFct, 
        double parameterRescalingProposal, 
        double hyperShape, double hyperScale,
        double initialAgreementParameter)
    {
      if (parameterRescalingProposal <= 1) throw new RuntimeException();
//      if (hyperParamMean <= 0) throw new RuntimeException();
      this.agreementParameter = initialAgreementParameter;
      this.suffStatExtractor = suffStatExtractor;
      this.lossFct = lossFct;
      this.a = parameterRescalingProposal;
      this.hyperShape = hyperShape;
      this.hyperScale = hyperScale;
    }
    public void update(DatasetType datasetBeingResampled, UnrootedTree newValue)
    {
      states.put(datasetBeingResampled, newValue);
      statesSS.put(datasetBeingResampled, suffStatExtractor.extract(newValue));
    }
    public void addToJointModel(
        final DatasetType datasetBeingResampled,
        final PhyloSampler ps)
    {
      if (Double.isNaN(this.temp)) 
        this.temp = ps.getTemperature();
      else if (this.temp != ps.getTemperature())
        throw new RuntimeException();
      final NonClockTreePrior basePrior = ps.getPrior();
      basePriors.put(datasetBeingResampled, basePrior);
//      baseSamplers.put(datasetBeingResampled, ps);
      final NonClockTreePrior conditionalPrior = new NonClockTreePrior() {
        private DatasetType other = null; 
        @Override public double logPriorDensity(UnrootedTree nct)
        {
          if (other == null) other = other(states.keySet(), datasetBeingResampled);
          return (basePriors.get(datasetBeingResampled).logPriorDensity(nct) + 
            logPriorAgreement(suffStatExtractor.extract(nct), statesSS.get(other))); ///temp; (NO! this one gets divided by PhyloSampler!)
        }
      };
      ps.setPrior(conditionalPrior);
    }
    public void sampleAgreementParameter(Random rand)
    {
      final double m = Sampling.nextDouble(rand, 1.0/a, a);
      final double newParam = m * agreementParameter;
      final double hyperDensityLogRatio = hyperLogDensity(newParam) - hyperLogDensity(agreementParameter);
      final double densityLogRatio = logJointPriorDensity(newParam) - logJointPriorDensity(agreementParameter);
      final double ratio = Math.min(1,Math.exp(Math.log(m) + (hyperDensityLogRatio + densityLogRatio)/temp));
      resamplingStats.addValue(ratio);
      if (rand.nextDouble() < ratio)
        agreementParameter = newParam;
    }
    public void increaseAgreementParam(double monotonicIncrease)
    {
      agreementParameter += monotonicIncrease;
    }
    private double hyperLogDensity(double x)
    {
      return (hyperShape - 1) * Math.log(x) - x / hyperScale;
      //return Sampling.exponentialDensity(hyperParamMean, x);
    }
    public double getCurrentDistance()
    {
      List<SuffStat> suffStats = CollUtils.list();
      for (DatasetType dt : basePriors.keySet())
        suffStats.add(suffStatExtractor.extract(states.get(dt)));
      return lossFct.loss(suffStats.get(0), suffStats.get(1));
    }
    private double logJointPriorDensity(double param)
    {
      List<SuffStat> suffStats = CollUtils.list();
      double sum = 0.0;
      for (DatasetType dt : basePriors.keySet())
      {
        sum += basePriors.get(dt).logPriorDensity(states.get(dt));
        suffStats.add(suffStatExtractor.extract(states.get(dt)));
      }
      return sum + logPriorAgreement(suffStats.get(0), suffStats.get(1), param);
    }
    private double logPriorAgreement(SuffStat ss1, SuffStat ss2)
    {
      return logPriorAgreement(ss1, ss2, agreementParameter);
    }
    private double logPriorAgreement(SuffStat ss1, SuffStat ss2, double param)
    {
      return - param * lossFct.loss(ss1, ss2);
    }
  } 
  public static interface NonClockTreeSuffStatExtractor<SuffStat>
  {
    public SuffStat extract(UnrootedTree nct);
  }
  public class CladeExtractor implements NonClockTreeSuffStatExtractor<Set<Set<Taxon>>>
  {
    @Override public Set<Set<Taxon>> extract(UnrootedTree nct) { return ml.filterClades(nct.clades()); }
  }
  
  public static DatasetType other(Collection<DatasetType> all, DatasetType one)
  {
    if (all.size() != 2) 
      throw new RuntimeException();
    for (DatasetType ds : all)
      if (!ds.equals(one))
        return ds;
    throw new RuntimeException();
  }
}
