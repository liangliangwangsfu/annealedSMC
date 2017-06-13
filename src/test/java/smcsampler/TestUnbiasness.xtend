package smcsampler

import org.junit.Test
import bayonet.distributions.ExhaustiveDebugRandom
import mains.AnnealedApproxFactory
import mains.HMM
import smcsampler.algo.AnnealedSMC
import smcsampler.algo.schedules.FixedTemperatureSchedule
import mains.HMM.SimpleTwoStates
import baselines.AnnealingTypeAlgorithm
import smcsampler.algo.AnnealedParticle
import org.junit.Assert
import baselines.SteppingStone

class TestUnbiasness {
  
  @Test
  def void testAIS() {
    testAnnealedSMC(0.0)
  }
  
  @Test
  def void testAnnealedSMC() {
    testAnnealedSMC(1.0)
  }
  
  def void testAnnealedSMC(double essThreshold) {
    val exhausiveRand = new ExhaustiveDebugRandom
    val algo = new AnnealedSMC => [
      random = exhausiveRand
      nSamplesPerTemperature = 2
      resamplingESSThreshold = essThreshold
      temperatureSchedule = new FixedTemperatureSchedule => [
        nTemperatures = 3
      ]
    ]
    checkExpectedZEstimate(algo, exhausiveRand) 
  }
  
//  @Test   // (known to be biased)
  def void testSteppingStone() {
    val exhausiveRand = new ExhaustiveDebugRandom
    val algo = new SteppingStone => [
      random = exhausiveRand
      nSamplesPerTemperature = 2
      nBurnInPerTemperature = 0
      temperatureSchedule = new FixedTemperatureSchedule => [
        nTemperatures = 2
      ]
    ]
    checkExpectedZEstimate(algo, exhausiveRand) 
  }
  
  def static void checkExpectedZEstimate(AnnealingTypeAlgorithm<AnnealedParticle> approximationAlgo, ExhaustiveDebugRandom exhausiveRand) {
    println(approximationAlgo.class.simpleName)
    val hmm = new HMM => [
      len = 2
      inPlace = false
      parameters = new SimpleTwoStates
    ]
    println("observations = " + hmm.observations)
    val factory = new AnnealedApproxFactory => [
      forbidOutputFiles = true
      model = hmm
      approximationAlgorithm = approximationAlgo
    ]
    var expectation = 0.0
    var nProgramTraces = 0
    var totalPr = 0.0
    while (exhausiveRand.hasNext) {
      val logZ = factory.buildAndRun.logNormEstimate
      expectation += Math.exp(logZ) * exhausiveRand.lastProbability
      totalPr += exhausiveRand.lastProbability
      nProgramTraces++
    }
    println("nProgramTraces = " + nProgramTraces)
    println("expectation = " + expectation)
    println("truth = " + Math.exp(hmm.exactLogZ))
    println("---")
    Assert.assertEquals(Math.exp(hmm.exactLogZ), expectation, 1e-10)
    
  }
}
