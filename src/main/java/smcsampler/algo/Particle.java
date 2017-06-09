package smcsampler.algo;

public interface Particle
{
  /**
   * pi_nextTemperature / pi_temperature (this)
   */
  double logDensityRatio(double temperature, double nextTemperature);
}