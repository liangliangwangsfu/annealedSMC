package smcsampler.algo;

import java.util.Arrays;
import java.util.Random;

import bayonet.smc.ParticlePopulation;
import bayonet.smc.ResamplingScheme;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import briefj.BriefParallel;
import smcsampler.algo.schedules.AdaptiveTemperatureSchedule;
import smcsampler.algo.schedules.TemperatureSchedule;

public class AnnealedSMC<P extends Particle>
{
  @Arg                    @DefaultValue("0.5")
  public double resamplingESSThreshold = 0.5;
  
  @Arg                                  @DefaultValue("AdaptiveTemperatureSchedule")
  public TemperatureSchedule temperatureSchedule = new AdaptiveTemperatureSchedule();
  
  @Arg                                         @DefaultValue("MULTINOMIAL")
  public ResamplingScheme resamplingScheme = ResamplingScheme.MULTINOMIAL;

  @Arg     @DefaultValue("1_000")
  public int nParticles = 1_000;

  @Arg   @DefaultValue("1")
  public int nThreads = 1;
  
  final Kernels<P> proposal;
  
  /**
   * @return The particle population at the last step
   */
  public ParticlePopulation<P> sample(Random random)
  {
    Random [] parallelRandomStreams = SMCStaticUtils.parallelRandomStreams(random, nParticles);
    ParticlePopulation<P> population = initialize(parallelRandomStreams);
    
    double temperature = 0.0;
    while (temperature < 1.0)
    {
      double nextTemperature = temperatureSchedule.nextTemperature(population, temperature); 
      population = propose(parallelRandomStreams, population, temperature, nextTemperature);
      if (population.getRelativeESS() < resamplingESSThreshold && nextTemperature < 1.0)
        population = resample(random, population); 
      temperature = nextTemperature;
//      System.out.println("T = " + temperature + ", ESS = " + population.getRelativeESS());
    }
    return population;
  }

  private ParticlePopulation<P> resample(Random random, ParticlePopulation<P> population)
  {
    population = population.resample(random, resamplingScheme);
    if (proposal.inPlace())
    {
      P previous = null;
      for (int i = 0; i < nParticles; i++)
      {
        P current = population.particles.get(i);
        if (current == previous)
          population.particles.set(i, proposal.deepCopy(current)); 
        previous = current;
      }
    }
    return population;
  }

  private ParticlePopulation<P> propose(Random [] randoms, final ParticlePopulation<P> currentPopulation, double temperature, double nextTemperature)
  {
    final boolean isInitial = currentPopulation == null;
    
    final double [] logWeights = new double[nParticles];
    @SuppressWarnings("unchecked")
    final P [] particles = (P[]) new Particle[nParticles];
    
    BriefParallel.process(nParticles, nThreads, particleIndex ->
    {
      P proposed = isInitial ?
        proposal.sampleInitial(randoms[particleIndex]) :
        proposal.sampleNext(randoms[particleIndex], currentPopulation.particles.get(particleIndex), nextTemperature);
      logWeights[particleIndex] = 
        (isInitial ? 0.0 : proposed.incrementalLogWeight(temperature, nextTemperature) + Math.log(currentPopulation.getNormalizedWeight(particleIndex)));
      particles[particleIndex] = proposed;
    });
    
    return ParticlePopulation.buildDestructivelyFromLogWeights(
        logWeights, 
        Arrays.asList(particles),
        isInitial ? 0.0 : currentPopulation.logScaling);
  }
  
  private ParticlePopulation<P> initialize(Random [] randoms)
  {
    return propose(randoms, null, Double.NaN, Double.NaN);
  }
  
  public AnnealedSMC(Kernels<P> proposal)
  {
    this.proposal = proposal;
  }
}
