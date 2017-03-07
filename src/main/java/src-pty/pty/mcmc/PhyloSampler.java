package pty.mcmc;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import nuts.io.IO;
import nuts.math.Sampling;
import nuts.util.CollUtils;
import nuts.util.EasyFormat;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import pty.UnrootedTree;
import pty.UnrootedTree.UnrootedTreeProcessor;
import pty.eval.Purity;
import pty.mcmc.ProposalDistribution.StochasticNearestNeighborInterchangeProposal;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import fig.exec.Execution;
import goblin.Taxon;

public class PhyloSampler
{
  // things that do not need to be init at each exec of sample
  private ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
  private PhyloSampler.Options phyloSamplerOptions = _defaultPhyloSamplerOptions;
  private PriorOptions priorOptions = _defaultPriorOptions;
  
  // things to initialize before each call of sample()
  private double temperature = Double.NaN;

  private boolean initialized = false;
  private UnrootedTreeState initialState = null;
  private UnrootedTreeState currentState = null;
  private NonClockTreePrior prior = null;
  private int iteration = 0;
  private SummaryStatistics globalRatioStatistics = null;
  private Map<String,SummaryStatistics> detailedRatioStat = null;
  private SummaryStatistics conditioningAcceptStat = null;
  private List<PhyloProcessor> processors = null;
  private RecordHighestLikelihood rhl = null;
  private LinkedList<ProposalDistribution> proposalDistributions = null;
  private boolean outputText = false;
  private String fileOutputPrefix = "";
  
  public void setOutputText(boolean v) { this.outputText = v; }

  public void init(UnrootedTreeState initialState)
  {
    this.iteration = 0;
    this.temperature = 1.0;
    this.initialState = initialState;
    this.currentState = initialState;
    this.prior = phyloSamplerOptions.prior.prior(priorOptions);
    this.globalRatioStatistics = new SummaryStatistics();
    this.conditioningAcceptStat = new SummaryStatistics();
    this.detailedRatioStat = new HashMap<String,SummaryStatistics>();
    this.proposalDistributions = new LinkedList<ProposalDistribution>();
    this.processors  = new ArrayList<PhyloProcessor>();
    this.rhl = new RecordHighestLikelihood();
    this.processors.add(rhl);
    initialized = true;
		// log(0); // TODO: uncomment this later. I commented this because it
		// caused problems when
		// nThread>1: multiple threads try to access the same file
		// 'samples-0.newick.gz'. --- L. Wang, Nov 19,
		// 2016
  }
  /**
   * Create a PhyloSampler with the same options, but different temperature
   * @param temperature
   * @return
   */
  public PhyloSampler createHeatedVersion(double temperature)
  {
    PhyloSampler ps = new PhyloSampler();
    ps.outputText = false;
    ps.proposalOptions = this.proposalOptions;
    ps.phyloSamplerOptions = this.phyloSamplerOptions;
    ps.priorOptions = this.priorOptions;
    ps.setConstraints(this.conditionedCluster); //    ps.conditionedClades = this.conditionedClades;
    // important: clone the tree/state 
    // o.w. can create multi-thread issues because of caching in NonClockTree's
    ps.init(this.initialState.deepClone()); 
    ps.temperature = temperature;
    return ps;
  }
  public UnrootedTreeState mle() { return rhl.argmax; }
  public double mleLogLikelihood() { return rhl.argmax != null ? mle().logLikelihood() : Double.NaN; }
  public double logLikelihood() { return currentState.logLikelihood(); }
  public void sampleManyTimes() // as many times as options.nIterations
  {
    if (outputText) LogInfo.track("Sampling " + phyloSamplerOptions.nIteration + " MCMC steps");
    while (iteration < phyloSamplerOptions.nIteration)
      sample();
    if (outputText) LogInfo.end_track();
  }
  // number of moves taken by this chain
  public int nIterations() { return iteration; }
  public void sample()
  {
    if (!initialized) throw new RuntimeException();
    sample(phyloSamplerOptions.rand);
    if ((iteration+1)%phyloSamplerOptions.logFrequency == 0)
      log(iteration);
    iteration++;
  }
  public double getMeanAcceptanceRatio() { return globalRatioStatistics.getMean(); }
  /**
   * Notes:
   * -- does not include the temperature division!
   * -- do not cache this!!
   */
  public double energy(UnrootedTreeState ncts) 
  { 
    final double value =  - ncts.logLikelihood() - prior.logPriorDensity(ncts.getNonClockTree());
    return value;
  }
  public double currentEnergy() { return energy(currentState); }
  public double energyOverTemperature(UnrootedTreeState ncts) 
  { 
    return (- ncts.logLikelihood() - prior.logPriorDensity(ncts.getNonClockTree()))/temperature;
  }
  public double currentEnergyOverTemperature() { return energyOverTemperature(currentState); }
  
  public static void sampleSwap(PhyloSampler ps1, PhyloSampler ps2, Random rand, SummaryStatistics ratioStats)
  {
    final double ratio = swapRatio(ps1, ps2);
    ratioStats.addValue(ratio);
    if (rand.nextDouble () < ratio)
      swapStates(ps1,ps2);
  }
  
  private static void swapStates(PhyloSampler ps1, PhyloSampler ps2)
  {
    final UnrootedTreeState s1 = ps1.currentState;
    final UnrootedTreeState s2 = ps2.currentState;
    ps1.currentState = s2;
    ps2.currentState = s1;
  }
  private static double swapRatio(PhyloSampler ps1, PhyloSampler ps2)
  {
    return Math.min(1.0, 
        Math.exp(
            (ps1.currentEnergy() - ps2.currentEnergy()) * 
            (1.0/ps1.getTemperature() - 1.0/ps2.getTemperature())));
  }
  private PrintWriter currentTreeOut = null;
  private int nTreesInCurrentFile = 0;
  private UnrootedTree last = null;
  private void log(int iter)
  {
    if (outputText) LogInfo.logs(this);
    if (temperature == 1.0) // only log samples if heated
      outputTree(currentState.getNonClockTree(), iter);
    // output to file sample, best sample, print statistics
    if (mle() != null && last != mle().getNonClockTree())
      outputMLETree(last = mle().getNonClockTree());
  }
  private void outputMLETree(UnrootedTree nct)
  {
    IO.writeToDisk(
        Execution.getFile(fileOutputPrefix + "mle.newick"), 
        nct.toNewick());
  }
  private void outputTree(UnrootedTree nct, int iter)
  {
    if (currentTreeOut == null || nTreesInCurrentFile >= phyloSamplerOptions.nTreesPerFile)
    {
      if (currentTreeOut != null) currentTreeOut.close();
      currentTreeOut = IOUtils.openOutHard(Execution.getFile(fileOutputPrefix + "samples-" + iter + ".newick.gz"));
      nTreesInCurrentFile = 0;
    }
    currentTreeOut.append(nct.toNewick() + "\n");
    nTreesInCurrentFile++;
  }
  public void closeFile() 
  {
    if (currentTreeOut != null)
    {
      currentTreeOut.close();
      currentTreeOut = null;
    }
  }
  
//  private Set<Set<Language>> conditionedClades = null;
//  private boolean gotInConditionedSet = false; No! not correct b/c of chain swaps!
  private Map<Taxon,String> conditionedCluster = null;
  
  public boolean isConditioning() { return conditionedCluster != null; } //conditionedClades != null; }
  
  public void  setConstraints(
      Map<Taxon,String> constraints)
  {
    this.conditionedCluster = constraints;
  }
//  {
//    constraints = new HashMap<Language,T>(constraints);
//    Set<Language> leaves = currentState.getNonClockTree().leavesSet();
//    constraints.keySet().retainAll(leaves);
//    
//    if (!constraints.keySet().equals(leaves))
//      throw new RuntimeException();
//    
//    
//    if (conditionedClades == null)
//      conditionedClades = new HashSet<Set<Language>>();
//    
//    conditionedClades.addAll(CollUtils.inducedPartition(constraints));
//    
//    checkConstraintsConsistent(conditionedClades);
//  }
  
//  private static void checkConstraintsConsistent(Set<Set<Language>> conditionedClades)
//  {
//    for (Set<Language> clade1 : conditionedClades)
//      for (Set<Language> clade2 : conditionedClades)
//        if (CollUtils.overlap(clade1, clade2))
//          throw new RuntimeException();
//  }
  
//  private double fractionOfConditionedCladesViolated(NonClockTree t)
//  {
//    final double condSize = conditionedClades.size();
//    return (condSize - CollUtils.inter(conditionedClades, t.clades()).size())/condSize;
//  }
  
  private double conditionedImpurity(UnrootedTree t)
  {
    return 1.0 - Purity.purity(t, conditionedCluster);
  }
  
  private boolean annealAccept(Random rand, UnrootedTreeState proposedState)
  {
    final boolean result = _annealAccept( rand,  proposedState);
    conditioningAcceptStat.addValue(result == true ? 1.0 : 0.0);
    return result;
  }
  
  public double getConditionalAnnealRatio() { return conditioningAcceptStat.getMean(); }
  public double getConditionalFraction() { return conditionedImpurity(currentState.getNonClockTree()); } //fractionOfConditionedCladesViolated(currentState.getNonClockTree()); }
  
  private boolean _annealAccept(Random rand, UnrootedTreeState proposedState)
  {
    if (!isConditioning()) return true;
    final double 
      newFractionOfConditionedCladesViolated = conditionedImpurity(proposedState.getNonClockTree());//fractionOfConditionedCladesViolated(proposedState.getNonClockTree());
    if (newFractionOfConditionedCladesViolated == 0.0) return true;
    final double 
      oldFractionOfConditionedCladesViolated = getConditionalFraction();
    if (oldFractionOfConditionedCladesViolated == 0.0) return false; //if the clades constraints were satisfied, dont break them
    // 
    if (newFractionOfConditionedCladesViolated <= oldFractionOfConditionedCladesViolated)
      return true;
    // still allow sometimes symm diff to increase a bit
    // conjecture: the space of trees with clade conditioning is irreducible with sNNI but
    // I found a counter example where it is not possible to get to it by NNI moves that monotonically increase fraction violated
    // so need this anneal-like procedure until you hit the conditioning set
    return rand.nextDouble() < phyloSamplerOptions.conditionAnneal;
  }
  
  
  public boolean nnisample(Random rand, UnorderedPair<Taxon,Taxon> selectedEdge)
  {
	boolean accept=false;	
    final double  multiplicativeBranchProposalScaling = 2.0;
    //MultiplicativeBranchProposal proposal=new MultiplicativeBranchProposal(multiplicativeBranchProposalScaling,false); 
    //proposal.selectedEdge = selectedEdge;     
    StochasticNearestNeighborInterchangeProposal proposal=new StochasticNearestNeighborInterchangeProposal(true, multiplicativeBranchProposalScaling); 
    proposal.selectedEdge = selectedEdge;    
    final Pair<UnrootedTree,Double> result = proposal.propose(currentState.getNonClockTree(), rand);
    if (result != null) // might happen e.g. when trying to do nni with 3 leaves
    {
      final UnrootedTreeState proposedState = currentState.copyAndChange(result.getFirst());
      final double logProposalRatio = result.getSecond();
      final double energyLogRatio = energy(proposedState) - currentEnergy();
      final double ratio = Math.min(1,Math.exp(logProposalRatio - energyLogRatio/temperature));
      if (Double.isNaN(ratio)) throw new RuntimeException();
//      globalRatioStatistics.addValue(ratio);
//      getDetailedRatioStatistics(proposal).addValue(ratio);
      //if (annealAccept(rand, proposedState) && rand.nextDouble() < ratio)
      if (rand.nextDouble() < ratio)
      {
        currentState = proposedState;
        accept=true;
//        System.out.println("accepted");
      }
      //for (PhyloProcessor processor : processors)
       // processor.process(currentState);      
    }
    return accept;
  }


  
  
  
  public boolean sample(Random rand)
  {
	boolean accept=false;
    final ProposalDistribution proposal = nextProposal(rand);
   
    final Pair<UnrootedTree,Double> result = proposal.propose(currentState.getNonClockTree(), rand);
    if (result != null) // might happen e.g. when trying to do nni with 3 leaves
    {
      final UnrootedTreeState proposedState = currentState.copyAndChange(result.getFirst());
      final double logProposalRatio = result.getSecond();
      final double energyLogRatio = energy(proposedState) - currentEnergy();
      final double ratio = Math.min(1,Math.exp(logProposalRatio - energyLogRatio/temperature));
      if (Double.isNaN(ratio)) throw new RuntimeException();
      globalRatioStatistics.addValue(ratio);
      getDetailedRatioStatistics(proposal).addValue(ratio);
      if (annealAccept(rand, proposedState) && rand.nextDouble() < ratio)
      {
        currentState = proposedState;
        accept=true;
      }      
      for (PhyloProcessor processor : processors)
        processor.process(currentState);      
    }
    return accept;
  }

  private ProposalDistribution nextProposal(Random rand)
  {
    if (proposalDistributions.isEmpty())
      proposalDistributions.addAll(ProposalDistribution.Util.proposalList(proposalOptions, initialState.getNonClockTree(), rand));
    return proposalDistributions.get(rand.nextInt(proposalDistributions.size()));
  }
  private SummaryStatistics getDetailedRatioStatistics(ProposalDistribution pd)
  {
    return CollUtils.getNoNull(detailedRatioStat, pd.description(), new SummaryStatistics());
  }
  public static interface PhyloProcessor
  {
    public void process(UnrootedTreeState ncts);
  }
  
  public static class PhyloProcessorAdaptor implements PhyloProcessor
  {
    public final UnrootedTreeProcessor rtp;
    public PhyloProcessorAdaptor(UnrootedTreeProcessor rtp)
    {
      this.rtp = rtp;
    }
    @Override
    public void process(UnrootedTreeState ncts)
    {
      rtp.process(ncts.getNonClockTree());
    }
  }
  
  public static class RecordHighestLikelihood implements PhyloProcessor
  {
    private UnrootedTreeState argmax = null;
    private double max = Double.NEGATIVE_INFINITY;
    public void process(UnrootedTreeState ncts)
    {
      if (ncts.logLikelihood() > max)
      {
        max = ncts.logLikelihood();
        argmax = ncts;
      }
    }
  }
  
  
  public static interface NonClockTreePrior
  {
    public double logPriorDensity(UnrootedTree nct);
  }
  
  public static class ImproperPrior implements NonClockTreePrior
  {
    @Override
    public double logPriorDensity(UnrootedTree nct) { return 0; }
  }
  
  public static class ExponentialPrior implements NonClockTreePrior
  {
    public final double meanParam;
    public ExponentialPrior(double meanParam)
    {
      this.meanParam = meanParam;
    }
    public double logPriorDensity(UnrootedTree nct)
    {
      double sum = 0.0;
      for (UnorderedPair<Taxon,Taxon> edge : nct.edges())
        sum += Sampling.exponentialLogDensity(meanParam, nct.branchLength(edge));
      return sum;
    }
  }
  @Override
  public String toString()
  {
    String result = toString("global", globalRatioStatistics) + " {\n";
    result += "\tBest log likelihood: " + mleLogLikelihood() + "\n";
    result += "\tCurrent log likelihood: " + currentState.logLikelihood() + "\n";
    for (Object key : detailedRatioStat.keySet())
      result += "\t" + toString(key, detailedRatioStat.get(key)) + "\n";
    return result + "}\n";
  }
  public String toString(Object descr, SummaryStatistics ss)
  {
    return "Number of " + descr + " sampling steps: " + ss.getN() + "\t"
      + "(Mean acceptance ratio: " + ss.getMean()+")";
  }
  public String detailedRatioToString()
  {
    StringBuilder result = new StringBuilder(" ");
    List<String> keys = new ArrayList<String>(detailedRatioStat.keySet());
    Collections.sort(keys);
    for (String key : keys)
      result.append(key + "=" + EasyFormat.fmt2(detailedRatioStat.get(key).getMean()) +
          " ");
    return result.toString();
  }

  public ProposalDistribution.Options getProposalOptions()
  {
    return proposalOptions;
  }

  public void setProposalOptions(ProposalDistribution.Options proposalOptions)
  {
    this.proposalOptions = proposalOptions;
  }

  public PhyloSampler.Options getPhyloSamplerOptions()
  {
    return phyloSamplerOptions;
  }

  public void setPhyloSamplerOptions(PhyloSampler.Options phyloSamplerOptions)
  {
    this.phyloSamplerOptions = phyloSamplerOptions;
  }

  public List<PhyloProcessor> getProcessors()
  {
    return processors;
  }

  public void setProcessors(List<PhyloProcessor> processors)
  {
    this.processors = processors;
  }
  
  public double getTemperature()
  {
    return temperature;
  }
  public void setTemperature(double temperature)
  {
    this.temperature = temperature;
  }
  public void setFileOutputPrefix(String fileOutputPrefix)
  {
    this.fileOutputPrefix = fileOutputPrefix;
  }
  public UnrootedTreeState getInitialState()
  {
    return initialState;
  }
  public UnrootedTreeState getCurrentState() { return currentState; }
  
  public static class Options
  {
    @Option public Prior prior = Prior.EXP;
    @Option public Random rand = new Random(1);
    @Option public int nIteration = 1000;
    @Option public int logFrequency = 100;
    @Option public int nTreesPerFile = 100;
    @Option public double conditionAnneal = 0.01;
  }
  public static class PriorOptions
  {
    @Option public double multiplicativeBranchFactor = 2.0;
  }
  public static PhyloSampler.Options _defaultPhyloSamplerOptions = new PhyloSampler.Options();
  public static PriorOptions _defaultPriorOptions = new PriorOptions();
  public static enum Prior
  {
    EXP {
      @Override
      public NonClockTreePrior prior(PriorOptions options)
      {
        //return new ExponentialPrior(options.multiplicativeBranchFactor);
    	  return new ExponentialPrior(10.0);
      }
    },
    IMPROPER {
      @Override
      public NonClockTreePrior prior(PriorOptions options)
      {
        return new ImproperPrior();
      }
    };
    abstract NonClockTreePrior prior(PriorOptions options);
  }
  public NonClockTreePrior getPrior()
  {
    return prior;
  }
  public void setPrior(NonClockTreePrior prior)
  {
    this.prior = prior;
  }


  

}
