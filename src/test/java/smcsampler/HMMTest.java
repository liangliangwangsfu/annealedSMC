package smcsampler;

import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import mains.HMM;
import mains.MeasureApproxFactory;
import smcsampler.algo.AnnealedSMC;
import smcsampler.algo.schedules.FixedTemperatureSchedule;

public class HMMTest
{
  @SuppressWarnings("unchecked")
  @Test
  public void testInPlace()
  {
    MeasureApproxFactory<?> factory = new MeasureApproxFactory<>();
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
  public void testConvergence()
  {
    MeasureApproxFactory<?> factory = new MeasureApproxFactory<>();
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
    
    Assert.assertEquals(exact, approx, Math.abs(exact) * 0.01);
  }
}
