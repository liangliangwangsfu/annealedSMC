package baselines;

import java.util.Random;

import bayonet.math.NumericalUtils;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import smcsampler.algo.Particle;
import smcsampler.algo.schedules.FixedTemperatureSchedule;
import smcsampler.algo.schedules.TemperatureSchedule;
import smcsampler.algo.Kernels;

public class SteppingStone<P extends Particle>
{
  @Arg              @DefaultValue("1_000")
  public int nMCMCPerTemperature = 1_000;
  
  @Arg                @DefaultValue("1_000")
  public int nBurnInPerTemperature = 1_000;
  
  @Arg               @DefaultValue("true")
  public boolean initFromPrevious = true;
  
  @Arg                                  @DefaultValue("FixedTemperatureSchedule")
  public TemperatureSchedule temperatureSchedule = new FixedTemperatureSchedule();
  
  
  final Kernels<P> proposal;
  
  public double estimateLogZ(Random random)
  {
    double sum = 0.0;
    double temperature = 0.0;
    P currentState = null;
    while (temperature < 1.0)
    {
      boolean isInit = temperature == 0.0;
      double nextTemperature = temperatureSchedule.nextTemperature(null, temperature); 
      double logSumAnnealedLikelihoods = Double.NEGATIVE_INFINITY;
      
      if (isInit)
        for (int i = 0; i < nBurnInPerTemperature; i++)
          currentState = proposal.sampleNext(random, currentState, temperature);
      else if (!initFromPrevious)
        currentState = proposal.sampleInitial(random);
      for (int i = 0; i < nMCMCPerTemperature; i++)
      {
        currentState = isInit ? 
          proposal.sampleInitial(random) :
          proposal.sampleNext(random, currentState, temperature);
        double logLikelihood = currentState.incrementalLogWeight(temperature, nextTemperature);
        logSumAnnealedLikelihoods = NumericalUtils.logAdd(logSumAnnealedLikelihoods, logLikelihood);
      }
      sum += logSumAnnealedLikelihoods - Math.log(nMCMCPerTemperature);
      temperature = nextTemperature;
    }
    return sum;
  }

  public SteppingStone(Kernels<P> proposal)
  {
    this.proposal = proposal;
  }
}
