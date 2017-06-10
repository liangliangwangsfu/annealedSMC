package smcsampler;

import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import baselines.SteppingStone;
import bayonet.smc.ParticlePopulation;
import mains.HMM;
import mains.AnnealedApproxFactory;
import smcsampler.algo.AnnealedSMC;
import smcsampler.algo.schedules.FixedTemperatureSchedule;

public class HMMTest
{
  @SuppressWarnings("unchecked")
  @Test
  public void testInPlace()
  {
    AnnealedApproxFactory<?> factory = new AnnealedApproxFactory<>();
    factory.forbidOutputFiles = true;
    @SuppressWarnings("rawtypes")
    AnnealedSMC algo = new AnnealedSMC<>();
    factory.approximationAlgorithm = algo;
    HMM hmm = new HMM();
    factory.model = hmm;
    hmm.inPlace = true;
    double v1 = factory.buildAndRun().logNormEstimate();
    hmm.inPlace = false;
    algo.random = new Random(1);
    double v2 = factory.buildAndRun().logNormEstimate();
    
    Assert.assertEquals(v1, v2, 0.0); 
  }
  
  @SuppressWarnings("unchecked")
  @Test
  public void testAnnealedConvergence()
  {
    AnnealedApproxFactory<?> factory = new AnnealedApproxFactory<>();
    factory.forbidOutputFiles = true;
    @SuppressWarnings("rawtypes")
    AnnealedSMC algo = new AnnealedSMC<>();
    algo.nParticles = 1000;
    factory.approximationAlgorithm = algo;
    HMM hmm = new HMM();
    factory.model = hmm;
    FixedTemperatureSchedule temp = new FixedTemperatureSchedule();
    algo.temperatureSchedule = temp;
    temp.nTemperatures = 1000;
    
    double approx = (factory.buildAndRun().logNormEstimate());
    double exact  = (hmm.getExactLogZ());
    System.out.println("Annealed SMC");
    System.out.println(approx);
    System.out.println(exact);
    Assert.assertEquals(exact, approx, Math.abs(exact) * 0.01);
  }
  
  @SuppressWarnings("unchecked")
  @Test
  public void testSteppingStoneConvergence()
  {
    AnnealedApproxFactory<?> factory = new AnnealedApproxFactory<>();
    factory.forbidOutputFiles = true;
    @SuppressWarnings("rawtypes")
    SteppingStone algo = new SteppingStone<>();
    algo.nMCMCPerTemperature = 1000;
    factory.approximationAlgorithm = algo;
    HMM hmm = new HMM();
    factory.model = hmm;
    FixedTemperatureSchedule temp = new FixedTemperatureSchedule();
    algo.temperatureSchedule = temp;
    temp.nTemperatures = 1000;
    
    ParticlePopulation<?> pop = factory.buildAndRun();
    System.out.println(pop.nParticles());
    double approx = (pop.logNormEstimate());
    double exact  = (hmm.getExactLogZ());
    System.out.println("Stepping Stone");
    System.out.println(approx);
    System.out.println(exact);
    Assert.assertEquals(exact, approx, Math.abs(exact) * 0.01);
  }
}
