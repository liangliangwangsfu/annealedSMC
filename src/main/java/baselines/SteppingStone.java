package baselines;

import java.util.Random;

import blang.inits.Arg;
import blang.inits.DefaultValue;
import smcsampler.algo.AdaptiveTemperatureSchedule;
import smcsampler.algo.FixedTemperatureSchedule;
import smcsampler.algo.Particle;
import smcsampler.algo.Proposal;
import smcsampler.algo.TemperatureSchedule;

public class SteppingStone<P extends Particle>
{
  @Arg              @DefaultValue("1_000")
  public int nMCMCPerTemperature = 1_000;
  
  @Arg                                  @DefaultValue("FixedTemperatureSchedule")
  public TemperatureSchedule temperatureSchedule = new FixedTemperatureSchedule();
  
  
  final Proposal<P> proposal;
  
  public double estimateLogZ(Random random)
  {
    double sum = 0.0;
    for (int i = 0; i < nTemperatures; i++)
    return sum;
  }

  public SteppingStone(Proposal<P> proposal)
  {
    this.proposal = proposal;
  }
}
