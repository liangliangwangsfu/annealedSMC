package mains;

import bayonet.smc.ParticlePopulation;
import blang.inits.Implementations;
import mains.AnnealedApproxFactory.EvaluationContext;
import smcsampler.algo.AnnealingKernels;
import smcsampler.algo.AnnealedParticle;

@Implementations({HMM.class})
interface Model<P extends AnnealedParticle>
{
  AnnealingKernels<P> kernels();
  
  void evaluate(ParticlePopulation<P> approximation, EvaluationContext context);
}