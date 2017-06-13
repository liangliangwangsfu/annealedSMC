package mains;

import java.util.Optional;
import java.util.Random;

import bayonet.smc.ParticlePopulation;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import mains.AnnealedApproxFactory.EvaluationContext;
import pty.UnrootedTree;
import smcsampler.algo.AnnealedLikelihoodParticle;
import smcsampler.algo.AnnealingKernels;

public class LegacyPhylo implements Model<AnnealedLikelihoodParticle<UnrootedTree>>
{
  @Arg    @DefaultValue("100")
  public    int nSites = 100;
  
  @Arg   @DefaultValue("10")
  public    int nTaxa = 10;
  
  @Arg               @DefaultValue("1")
  public Random random = new Random(1);

  @Override
  public AnnealingKernels<AnnealedLikelihoodParticle<UnrootedTree>> kernels()
  {
    // TODO 
    throw new RuntimeException();
  }

  @Override
  public void evaluate(ParticlePopulation<AnnealedLikelihoodParticle<UnrootedTree>> approximation,
      EvaluationContext context)
  {
    double approx = approximation.logNormEstimate();
    context.registerLogZ(approx, Optional.empty());
  }

}
