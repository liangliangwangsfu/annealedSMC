package smcsampler;


import org.junit.Assert;
import org.junit.Test;

import baselines.AnnealingTypeAlgorithm;
import baselines.SteppingStone;
import bayonet.distributions.Random;
import mains.HMM;
import mains.AnnealedApproxFactory;
import smcsampler.algo.AnnealedSMC;
import smcsampler.algo.schedules.FixedTemperatureSchedule;

public class TestConvergence
{

  @SuppressWarnings({ "rawtypes", "unchecked" })
  public static void test(AnnealingTypeAlgorithm<?> algo)
  {
    AnnealedApproxFactory factory = new AnnealedApproxFactory<>();
    factory.forbidOutputFiles = true;
    HMM hmm = new HMM();
    hmm.parameters = new HMM.SimpleTwoStates();
    hmm.len = 3;
    factory.model = hmm;
    factory.approximationAlgorithm = algo;
    
    double approx = (factory.buildAndRun().logNormEstimate());
    double exact  = (hmm.getExactLogZ());
    System.out.println(algo.getClass().getSimpleName());
    System.out.println(approx);
    System.out.println(exact);
    System.out.println("---");
    Assert.assertEquals(exact, approx, Math.abs(exact) * 0.001);
  }
  
  @SuppressWarnings("unchecked")
  @Test
  public void testAnnealedConvergence()
  {
    @SuppressWarnings("rawtypes")
    AnnealedSMC algo = new AnnealedSMC<>();
    algo.nSamplesPerTemperature = 1000;
    FixedTemperatureSchedule temp = new FixedTemperatureSchedule();
    algo.temperatureSchedule = temp;
    temp.nTemperatures = 1000;
    test(algo);
  }
  
  @SuppressWarnings("unchecked")
  @Test
  public void testSteppingStoneConvergence()
  {
    @SuppressWarnings("rawtypes")
    SteppingStone algo = new SteppingStone<>();
    algo.nSamplesPerTemperature = 1000;
    FixedTemperatureSchedule temp = new FixedTemperatureSchedule();
    algo.temperatureSchedule = temp;
    temp.nTemperatures = 1000;
    test(algo);
  }
  
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
  
}
