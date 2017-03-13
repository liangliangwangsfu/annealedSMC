package smcsampler;
import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import conifer.Phylogeny;
//import conifer.msa.TreeMSAState;
import conifer.particle.PhyloParticle;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.RootedTree.RootedTreeProcessor;
import pty.UnrootedTree.UnrootedTreeProcessor;
import pty.smc.LazyPCS;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter.ParticleProcessor;

import ev.ex.NJPState;
//import ev.ex.NJPState;
//import ev.ex.NJStateKernel;
import ev.poi.PoissonModel;
import ev.poi.PoissonSampleProcessor;
import ev.poi.SampleContext;
import ev.to.NJ;
import fig.basic.LogInfo;
import fig.basic.UnorderedPair;
import goblin.Taxon;
import nuts.util.Counter;


public class TreeDistancesProcessor 
  implements 
    PoissonSampleProcessor,
    RootedTreeProcessor,
    UnrootedTreeProcessor,
    ParticleProcessor,
    Serializable, ParticleFilterSMCSampler.ParticleProcessor
{
  private static final long serialVersionUID = 1L;
  private final Counter<UnorderedPair<Taxon,Taxon>> _meanDistances = new Counter<UnorderedPair<Taxon,Taxon>>();
  private double norm = 0.0;
  private List<Taxon> taxa = null;
  private double bestLogLikelihood = Double.NEGATIVE_INFINITY;
  private Phylogeny best = null; //PhyloParticle best = null;
  
  public double lastLogLikelihood;
  public Phylogeny lastPhylo;
  
  public Phylogeny getMode()
  {
    return best;
  }
  
  public double getBestLogLikelihood() { return bestLogLikelihood; }

  public Counter<UnorderedPair<Taxon,Taxon>> getUnrootedCladesPosterior()
  {
    Counter<UnorderedPair<Taxon,Taxon>> result = new Counter<UnorderedPair<Taxon,Taxon>>();
    for (UnorderedPair<Taxon,Taxon> key : _meanDistances.keySet())
      result.setCount(key, _meanDistances.getCount(key) / norm);
    return result;
  }
  
  public UnrootedTree getConsensus() { return getConsensus(false); }
  public UnrootedTree getConsensus(boolean useNinja)
  {
 //   UnrootedTree njTree = useNinja ?
//        ev.to.Ninja.inferTreeFromDistances(getUnrootedCladesPosterior()) :
//        NJ.inferTree(getUnrootedCladesPosterior());
        UnrootedTree njTree = NJ.inferTree(getUnrootedCladesPosterior());
    return UnrootedTree.removeZeroes(njTree);
  }
  @Override
  public void process(PoissonModel sample, SampleContext context)
  {
    process(sample.currentUnrooted());
  }
  @Override
  public void process(RootedTree rt)
  {
    process(UnrootedTree.fromRooted(rt));
  }
  @Override
  public void process(UnrootedTree ut)
  {
    norm += 1.0;
    _process(ut, 1.0);
    
  }
  @Override
  public void process(Object state, double weight)
  {
    if (state instanceof LazyPCS)
      state = ((LazyPCS) state).getState();
    if (state instanceof PartialCoalescentState)
    {
      process(UnrootedTree.fromRooted(((PartialCoalescentState)state).getFullCoalescentState()), weight);
      double currentLogLL = ((PartialCoalescentState) state).logLikelihood();
      if (currentLogLL > bestLogLikelihood)
      {
        best = ((PartialCoalescentState) state).getFullCoalescentState();
        bestLogLikelihood = currentLogLL;
      }
    }
    else if (state instanceof PhyloParticle)
    {
//      the best likelihood F1 does not add up if running:
//      -phy.data.genEvol.model PIP  -phy.data.genEvol.pipParamType  LEN_INT -phy.infEvo.model  PIP -phy.infEvo.pipParamType   LEN_INT -phy.infEvo.pipParam1  20 -phy.infEvo.pipParam2  .1  -phy.exp.experiments PARAM_MONITOR MSA_HELDOUT TREE_HELDOUT -phy.mc.algorithm TEMPERED_MCMC -phy.prop.proposalModels  INFORMED_LOCAL_MSA -phy.prop.weights 1.0 -phy.data.genTree.nTaxa 20 -phy.prop.nGenerations MAX -phy.masterRand.masterRandom 991 -dataFile data/gutell/raw/5S.3.alnfasta -dataModel GENERATED -phy.prop.forceTreeMSAParticles true -generationsPerProcessPeriod 10 -processesPerMonitorPeriod 100 -phy.data.genEvol.tkfParam2 0.5 -phy.data.genTree.rate 0.5 -phy.prop.treeInitMethod  DATA  -phy.data.genEvol.tkfParam1 1000 -phy.prop.msaInitMethod DATA -phy.prop.initmsafromdataevenwhen true -phy.data.genEvol.pipParam1 20 -phy.data.genEvol.pipParam2 .1 -phy.prop.initTreeFromDataEvenWhenHeldout -phy.infTree.rate 0.5
      
      
      PhyloParticle pp = (PhyloParticle) state;
      process(pp.getPhylogeny().getUnrooted(), weight);
      double currentLogLL = pp.getLogLikelihood();
      if (currentLogLL > bestLogLikelihood)
      {
        LogInfo.logsForce("better LL found: " + bestLogLikelihood + " -> [ ll(" + currentLogLL + ") + lp(" + pp.getLogPrior() + ") = " + (currentLogLL + pp.getLogPrior()) + " ]");
        LogInfo.logsForce("ll" + currentLogLL);
//        LogInfo.logsForce("n edges = " + ((TreeMSAState) pp).msa.edges().size());
//        LogInfo.logsForce("hash = " + set(((TreeMSAState) pp).msa.edges()).hashCode());
//        LogInfo.logsForce("f1 = " + MSAPoset.edgeF1(ParallelTemperedMCMC.getRef(), ((TreeMSAState) pp).msa));
        best = pp.getPhylogeny();
        bestLogLikelihood = currentLogLL;
      }
      
      lastLogLikelihood = currentLogLL;
      lastPhylo = pp.getPhylogeny();
    }
    else if (state instanceof NJPState)
      process(((NJPState) state).pcs, weight);
//    else 
//    	if (state instanceof NJPState2)
//    	      process(((NJPState2) state).pcs, weight);
    else
      throw new RuntimeException();
//    norm += weight;
//    _process(UnrootedTree.fromRooted(state.getFullCoalescentState()), weight);
  }
  public void process(UnrootedTree t, double w)
  {
    norm += w;
    _process(t, w);
  }
  private void _process(UnrootedTree ut, double w)
  {
    if (taxa == null)
      taxa = ut.leaves();
    
    
//    long beg = System.currentTimeMillis();
//    System.out.println("Test---remove me");
//    ut.allTotalBranchLengthDistances();
//    System.out.println(System.currentTimeMillis() - beg);
    // end remove me
    
    Counter<UnorderedPair<Taxon,Taxon>> c = ut.allTotalBranchLengthDistances();
    
    synchronized (this)
    {
      for (Entry<UnorderedPair<Taxon,Taxon>,Double> keyPair : c.entries.entrySet())
      {
        double count = keyPair.getValue();
        UnorderedPair<Taxon, Taxon> key =keyPair.getKey();
        _meanDistances.incrementCount(key, w * count);
      }
    }
//    _meanDistances.incrementAll();
    
//    for (int i = 0 ; i < taxa.size(); i++)
//      for (int j = i+1; j< taxa.size() ;j ++)
//      {
//        Taxon t1 = taxa.get(i), t2 = taxa.get(j);
//        _meanDistances.incrementCount(new UnorderedPair<Taxon,Taxon>(t1,t2), ut.totalBranchLengthDistance(t1,t2) * w);
//      }
//    System.out.println(System.currentTimeMillis() - beg);
  }

  public Counter<UnorderedPair<Taxon, Taxon>> getMeanDistances()
  {
    return new Counter<UnorderedPair<Taxon,Taxon>>(_meanDistances);
//    return null;
  }
}