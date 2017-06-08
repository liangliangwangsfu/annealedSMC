package smcsampler;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jgrapht.UndirectedGraph;
import org.junit.Test;

import bayonet.distributions.Bernoulli;
import bayonet.distributions.Multinomial;
import bayonet.graphs.GraphUtils;
import bayonet.marginal.DiscreteFactorGraph;
import bayonet.marginal.algo.SumProduct;
import smcsampler.algo.AnnealedLikelihoodParticle;
import smcsampler.algo.AnnealedSMC;
import smcsampler.algo.Kernels;
import smcsampler.algo.schedules.FixedTemperatureSchedule;

public class HMMTest
{
  boolean inPlace = false;
  int len = 3;
  
  double [][] transitionPrs = new double[][]{{0.8,0.15,0.05},{0.02,0.93,0.05},{0.15,0.15,0.7}};
  int nLatents = transitionPrs.length;
  double [][] emissionPrs = transitionPrs;
  double [] initialPrs = new double[]{0.25, 0.25, 0.5};

  Random random = new Random(1);
  List<Integer> observations = generateData(random);
  
  // compute true logZ using discrete HMM
  double trueLogZ;
  {
    DiscreteFactorGraph<Integer> dfg = createHMM(observations);
    SumProduct<Integer> sp = new SumProduct<>(dfg);
    trueLogZ = sp.logNormalization();
  }
  
  @Test
  public void testInPlace()
  {
    AnnealedSMC<AnnealedLikelihoodParticle<List<Integer>>> smc = 
        new AnnealedSMC<>(new HMMProposal());
    smc.nParticles = 100;
    
    inPlace = false;
    System.out.println("notInPlace = " + smc.sample(new Random(1)).logNormEstimate());
    inPlace = true;
    System.out.println("inPlace    = " + smc.sample(new Random(1)).logNormEstimate());
  }
  
//  @Test
//  public void test()
//  {
//    AnnealedSMC<AnnealedLikelihoodParticle<List<Integer>>> smc = 
//        new AnnealedSMC<>(new HMMProposal());
//    smc.nParticles = 10_000;
//    FixedTemperatureSchedule schedule = new FixedTemperatureSchedule();
//    schedule.nTemperatures = 10000;
//    smc.temperatureSchedule = schedule;
//    
//    System.out.println("smc   = " + smc.sample(random).logNormEstimate());
//    System.out.println("truth = " + trueLogZ);
//  }
  
  double logLikelihood(List<Integer> latents)
  {
    double sum = 0.0;
    for (int i = 0; i < len; i++)
      sum += Math.log(emissionPrs[latents.get(i)][observations.get(i)]);
    return sum;
  }
  
  double logPrior(List<Integer> latents)
  {
    double sum = Math.log(initialPrs[latents.get(0)]);
    for (int i = 1 ; i < len; i++)
      sum += Math.log(transitionPrs[latents.get(i-1)][latents.get(i)]);
    return sum;
  }
  
  class HMMProposal implements Kernels<AnnealedLikelihoodParticle<List<Integer>>>
  {
    @Override
    public AnnealedLikelihoodParticle<List<Integer>> sampleInitial(Random random)
    {
      List<Integer> fromPrior = samplePrior(random);
      return new AnnealedLikelihoodParticle<List<Integer>>(logLikelihood(fromPrior), fromPrior);
    }

    @Override
    public AnnealedLikelihoodParticle<List<Integer>> sampleNext(Random random,
        AnnealedLikelihoodParticle<List<Integer>> current, double temperature)
    {
      // pick one state at random
      int index = random.nextInt(len);
      int oldState = current.contents.get(index);
      List<Integer> states = inPlace ? current.contents : new ArrayList<>(current.contents);
      double logDenom = logPrior(states) + temperature * current.logLikelihood;
      
      int newState = random.nextInt(nLatents);
      states.set(index, newState);
      double logNum   = logPrior(states) + temperature * logLikelihood(states);
       
      double acceptPr = Math.min(1.0, Math.exp(logNum - logDenom));
      if (Bernoulli.generate(random, acceptPr))
        ;
      else
        states.set(index, oldState);
      
      return new AnnealedLikelihoodParticle<List<Integer>>(logLikelihood(states), states);
    }

    @Override
    public boolean inPlace()
    {
      return inPlace;
    }

    @Override
    public AnnealedLikelihoodParticle<List<Integer>> deepCopy(AnnealedLikelihoodParticle<List<Integer>> particle)
    {
      return new AnnealedLikelihoodParticle<>(particle.logLikelihood, new ArrayList<>(particle.contents));
    }
  }
  
  List<Integer> samplePrior(Random random)
  {
    List<Integer> prior = new ArrayList<>();
    int latent = Multinomial.sampleMultinomial(random, initialPrs);
    prior.add(latent);
    for (int i = 1; i < len; i++)
    {
      latent = Multinomial.sampleMultinomial(random, transitionPrs[latent]);
      prior.add(latent);
    }
    return prior;
  }
  
  List<Integer> generateData(Random random)
  {
    List<Integer> prior = samplePrior(random);
    List<Integer> result = new ArrayList<>();
    for (int i = 0; i < len; i++)
    {
      int latent = prior.get(i);
      result.add(Multinomial.sampleMultinomial(random, emissionPrs[latent]));
    }
    return result;
  }
  
  DiscreteFactorGraph<Integer> createHMM(List<Integer> observations)
  {
    UndirectedGraph<Integer, ?> topology = GraphUtils.createChainTopology(len);
    DiscreteFactorGraph<Integer> result = new DiscreteFactorGraph<Integer>(topology);
    
    // initial distribution
    result.setUnary(0, new double[][]{initialPrs});
    
    // transition
    for (int i = 0; i < len-1; i++)
      result.setBinary(i, i+1, transitionPrs);
    
    // observations
    for (int i = 0; i < len; i++)
    {
      int currentObs = observations.get(i);
      double [] curEmissionPrs = new double[initialPrs.length];
      for (int s = 0; s < initialPrs.length; s++)
        curEmissionPrs[s] = emissionPrs[s][currentObs];
      result.unaryTimesEqual(i, new double[][]{curEmissionPrs});
    }
    
    return result;
  }
}
