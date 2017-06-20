package baselines;

import java.util.ArrayList;
import java.util.List;

import bayonet.distributions.Random;
import bayonet.math.NumericalUtils;
import bayonet.smc.ParticlePopulation;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import smcsampler.algo.AnnealedParticle;
import smcsampler.algo.schedules.FixedTemperatureSchedule;
import smcsampler.algo.schedules.TemperatureSchedule;
import smcsampler.algo.AnnealingKernels;

public class SteppingStone<P extends AnnealedParticle> implements AnnealingTypeAlgorithm<P>
{
  @Arg                 @DefaultValue("1_000")
  public int nSamplesPerTemperature = 1_000;
  
  @Arg                @DefaultValue("1_000")
  public int nBurnInPerTemperature = 1_000;
  
  @Arg               @DefaultValue("true")
  public boolean initFromPrevious = true;
  
  @Arg                                  @DefaultValue("FixedTemperatureSchedule")
  public TemperatureSchedule temperatureSchedule = new FixedTemperatureSchedule();
  
  @Arg               @DefaultValue("1")
  public Random random = new Random(1);
  
  AnnealingKernels<P> kernels;
  
  @Override
  public ParticlePopulation<P> getApproximation()
  {
    double sum = 0.0;
    double temperature = 0.0;
    P currentState = kernels.sampleInitial(random);
    while (temperature < 1.0)
    {
      boolean isInit = temperature == 0.0;
      double nextTemperature = temperatureSchedule.nextTemperature(null, temperature); 
      double logSumAnnealedLikelihoods = Double.NEGATIVE_INFINITY;
      
      if (!initFromPrevious)
        currentState = kernels.sampleInitial(random);
      
      for (int i = 0; i < nBurnInPerTemperature; i++)
        currentState = kernels.sampleNext(random, currentState, temperature);
      
      for (int i = 0; i < nSamplesPerTemperature; i++)
      {
        currentState = isInit ? 
          kernels.sampleInitial(random) :
          kernels.sampleNext(random, currentState, temperature);
        double logLikelihood = currentState.logDensityRatio(temperature, nextTemperature);
        logSumAnnealedLikelihoods = NumericalUtils.logAdd(logSumAnnealedLikelihoods, logLikelihood);
      }
      sum += logSumAnnealedLikelihoods - Math.log(nSamplesPerTemperature);
      temperature = nextTemperature;
    }
    List<P> particles = new ArrayList<>();
    for (int i = 0; i < nSamplesPerTemperature; i++)
    {
      currentState = kernels.sampleNext(random, currentState, temperature);
      particles.add(kernels.inPlace() ? kernels.deepCopy(currentState) : currentState);
    }
    return ParticlePopulation.buildEquallyWeighted(particles, sum + Math.log(nSamplesPerTemperature)); // + log(n) b/c logZ estimate add - log(n) in ParticlePopulation
  }

  @Override
  public void setKernels(AnnealingKernels<P> kernels)
  {
    this.kernels = kernels;
  }
}
