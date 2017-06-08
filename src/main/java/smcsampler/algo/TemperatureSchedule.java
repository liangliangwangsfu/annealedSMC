package smcsampler.algo;

import bayonet.smc.ParticlePopulation;
import blang.inits.Implementations;

@Implementations({AdaptiveTemperatureSchedule.class, FixedTemperatureSchedule.class}) 
public interface TemperatureSchedule
{
  double nextTemperature(ParticlePopulation<? extends Particle> population, double temperature);
}