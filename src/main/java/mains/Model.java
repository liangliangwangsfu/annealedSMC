package mains;

import bayonet.smc.ParticlePopulation;
import blang.inits.Implementations;
import mains.MeasureApproxFactory.EvaluationContext;
import smcsampler.algo.Kernels;
import smcsampler.algo.Particle;

@Implementations({HMM.class})
interface Model<P extends Particle>
{
  Kernels<P> kernels();
  
  void evaluate(ParticlePopulation<P> approximation, EvaluationContext context);
}