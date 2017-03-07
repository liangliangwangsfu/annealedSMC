package pty.smc;
import java.util.*;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.eval.SymmetricDiff;
import pty.io.TreeEvaluator;
import pty.smc.ParticleFilter.StoreProcessor;

import ev.poi.processors.TreeDistancesProcessor;
import gep.util.OutputManager;
import goblin.Taxon;
import nuts.math.Sampling;
import static nuts.util.CollUtils.*;


public class PMCMC
{
  private final ParticleFilter<PartialCoalescentState> pf;
  private final NCPriorPriorCondKernel kernel;
  private final TreeDistancesProcessor tdp;
  private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
  private RootedTree currentSample = null;
  
  public static OutputManager outMan = new OutputManager();
  
  public PMCMC(ParticleFilter<PartialCoalescentState> pf,
      NCPriorPriorCondKernel kernel, TreeDistancesProcessor tdp)
  {
    this.pf = pf;
    this.kernel = kernel;
    this.tdp = tdp;
  }

  public void next(Random rand)
  {
    RootedTree previousSample = currentSample;
    // log the sparsity of current mask
    final double currentSparsity = (currentSample == null ? 0.0 : ((double)kernel.getNMasks()) / ((double)currentSample.topology().nLeaves()-3));
    // sample from PF
    StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();
    double acceptPr = 0.0;
    try 
    {
      pf.sample(kernel, pro);
      // compute accept/reject
  
      final double logRatio = pf.estimateNormalizer() - previousLogLLEstimate;
  
      acceptPr = Math.min(1, Math.exp(logRatio));
      if (currentSample != null &&  Double.isInfinite(acceptPr))
        throw new RuntimeException();
      final boolean accept = Sampling.sampleBern(acceptPr, rand);
      if (accept)
      {
        // sample from sample
        PartialCoalescentState sampled = pro.sample(rand);
        currentSample = sampled.getFullCoalescentState();
        // set ll!
        previousLogLLEstimate = pf.estimateNormalizer();
      }
    }
    catch (Exception e)
    {
      // total fail!
    }
    // update tdp
    tdp.process(currentSample);
    // log some stats
    final int tSize = currentSample.topology().nLeaves();
    outMan.write("PMCMC", 
        "treeSize", tSize, 
        "acceptPr", acceptPr, 
        "maskSparsity", currentSparsity,
        "rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)));
    // sample mask from sample from sample
    Set<Set<Taxon>> newMask = set();
    final double currentRetProb = rand.nextDouble();
    // iterate, ignore trivials, sample berns
    Set<Set<Taxon>> clades = UnrootedTree.fromRooted(currentSample).unRootedClades();
    Set<Taxon> allTaxa = SymmetricDiff.allLeaves(clades);
    Set<Set<Taxon>> simplified = set(clades);
    for (Set<Taxon> clade : clades)
    {
      Set<Taxon> complement = set(allTaxa);
      complement.removeAll(clade);
      if (simplified.contains(complement) && simplified.contains(clade))
        simplified.remove(complement);
    }
    if (simplified.size() != clades.size() /2 )
      throw new RuntimeException();
    for (Set<Taxon> clade : simplified)
      if (nonTrivial(clade, allTaxa.size()))
        if (Sampling.sampleBern(currentRetProb, rand))
          newMask.add(clade);
    // set it (for next time)
    kernel.setConditioning(newMask);
  }

  private boolean nonTrivial(Set<Taxon> clade, int nLeaves)
  {
    if (clade.size() <= 1 || clade.size() >= nLeaves - 1) return false;
    else return true;
  }
}
