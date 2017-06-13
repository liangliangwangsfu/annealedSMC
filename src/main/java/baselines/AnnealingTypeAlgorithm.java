package baselines;

import bayonet.smc.ParticlePopulation;
import blang.inits.Implementations;
import smcsampler.algo.AnnealedSMC;
import smcsampler.algo.AnnealingKernels;
import smcsampler.algo.AnnealedParticle;

@Implementations({AnnealedSMC.class, SteppingStone.class})
public interface AnnealingTypeAlgorithm<P extends AnnealedParticle>
{
  void setKernels(AnnealingKernels<P> kernels);
  ParticlePopulation<P> getApproximation();
}
