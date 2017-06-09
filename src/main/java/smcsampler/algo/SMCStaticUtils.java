package smcsampler.algo;

import java.util.Random;
import java.util.SplittableRandom;

import bayonet.distributions.Multinomial;
import bayonet.smc.ParticlePopulation;

public class SMCStaticUtils
{
  public static double relativeESS(ParticlePopulation<? extends Particle> population, double temperature, double nextTemperature, boolean conditional)
  {
    double [] incrementalWeights = incrementalWeights(population, temperature, nextTemperature);
    double 
      numerator   = 0.0,
      denominator = 0.0;
    for (int i = 0; i < population.nParticles(); i++)
    {
      double factor = population.getNormalizedWeight(i) * incrementalWeights[i];
      numerator   += factor;
      denominator += factor * (conditional ? incrementalWeights[i] : factor);
    }
    return numerator * numerator / denominator;
  }
  
  public static double[] incrementalWeights(ParticlePopulation<? extends Particle> population, double temperature,
      double nextTemperature)
  {
    double [] result = new double[population.nParticles()];
    for (int i = 0; i < population.nParticles(); i++)
      result[i] = population.particles.get(i).logDensityRatio(temperature, nextTemperature);
    Multinomial.expNormalize(result);
    return result;
  }
  
  public static Random [] parallelRandomStreams(Random seed, int nStreams)
  {
    Random[] result = new Random[nStreams]; 
    SplittableRandom splitRandom = new SplittableRandom(seed.nextLong());
    for (int i = 0; i < nStreams; i++)
      result[i] = new Random(splitRandom.split().nextLong());
    return result;
  }
}
