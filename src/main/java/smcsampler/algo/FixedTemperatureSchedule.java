package smcsampler.algo;

import bayonet.smc.ParticlePopulation;
import blang.inits.Arg;
import blang.inits.DefaultValue;

public class FixedTemperatureSchedule implements TemperatureSchedule
{
  @Arg        @DefaultValue("100")
  public int nTemperatures = 100;

  @Override
  public double nextTemperature(ParticlePopulation<? extends Particle> population, double temperature)
  {
    return Math.min(1.0, temperature + 1.0 / ((double) nTemperatures));
  }
}