package smcsampler.algo.schedules;

import bayonet.smc.ParticlePopulation;
import blang.inits.Implementations;
import smcsampler.algo.AnnealedParticle;

@Implementations({AdaptiveTemperatureSchedule.class, FixedTemperatureSchedule.class}) 
public interface TemperatureSchedule
{
  double nextTemperature(ParticlePopulation<? extends AnnealedParticle> population, double temperature);
}