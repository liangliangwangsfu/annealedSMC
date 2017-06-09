package baselines;

import bayonet.smc.ParticlePopulation;
import blang.inits.Implementations;
import smcsampler.algo.AnnealedSMC;
import smcsampler.algo.Kernels;
import smcsampler.algo.Particle;

@Implementations({AnnealedSMC.class})
public interface MeasureApproximation<P extends Particle>
{
  void setKernels(Kernels<P> kernels);
  ParticlePopulation<P> getApproximation();
}
