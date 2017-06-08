package smcsampler.algo;

public interface Particle
{
  double incrementalLogWeight(double temperature, double nextTemperature);
}