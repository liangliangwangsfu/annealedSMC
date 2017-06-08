package smcsampler.algo;

/**
 * Sequences of distribution of the form prior x like^temperature
 */
public class AnnealedLikelihoodParticle<P> implements Particle
{
  public final double logLikelihood;
  public final P contents;
  
  public AnnealedLikelihoodParticle(double logLikelihood, P contents)
  {
    this.logLikelihood = logLikelihood;
    this.contents = contents;
  }

  @Override
  public double incrementalLogWeight(double temperature, double nextTemperature)
  {
    return (nextTemperature - temperature) * logLikelihood;
  }
}