package mains;
 
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.jgrapht.UndirectedGraph;

import bayonet.distributions.Multinomial;
import bayonet.distributions.Random;
import bayonet.graphs.GraphUtils;
import bayonet.marginal.DiscreteFactorGraph;
import bayonet.marginal.algo.SumProduct;
import bayonet.smc.ParticlePopulation;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import blang.inits.DesignatedConstructor;
import blang.inits.Implementations;
import mains.AnnealedApproxFactory.EvaluationContext;
import smcsampler.algo.AnnealedLikelihoodParticle;
import smcsampler.algo.AnnealingKernels;

public class HMM implements Model<AnnealedLikelihoodParticle<List<Integer>>>
{
  @Arg                @DefaultValue("SimpleThreeStates")
  public Parameters parameters = new SimpleThreeStates();
  
  @Arg      @DefaultValue("false")
  public boolean inPlace = false;
  
  @Arg @DefaultValue("10")
  public    int len = 10;
  
  @Arg               @DefaultValue("1")
  public Random random = new Random(1);
  
  @Arg                       @DefaultValue("false")
  public boolean setAllObservationsToZero = false;
  
  @Override
  public AnnealingKernels<AnnealedLikelihoodParticle<List<Integer>>> kernels()
  {
    initCache();
    return new HMMKernels();
  }

  @Override
  public void evaluate(ParticlePopulation<AnnealedLikelihoodParticle<List<Integer>>> approximation,
      EvaluationContext context)
  {
    initCache();
    double approx = approximation.logNormEstimate();
    context.registerLogZ(approx, Optional.of(exactLogZ));
  }
  
  public List<Integer> observations()
  {
    initCache();
    return observations;
  }
  
  void initCache()
  {
    if (initialized)
      return;
    transitionPrs = parameters.transitionPrs();
    nLatents = transitionPrs.length;
    emissionPrs = parameters.emissionPrs();
    initialPrs = parameters.initialPrs();
    observations = generateData(random);
    DiscreteFactorGraph<Integer> dfg = createHMM(observations);
    SumProduct<Integer> sp = new SumProduct<>(dfg);
    exactLogZ = sp.logNormalization();
    
    initialized = true;
  }
  
  boolean initialized = false;
  double [][] transitionPrs = null;
  int nLatents = -1;
  double [][] emissionPrs = null;
  double [] initialPrs = null;
  List<Integer> observations = null;
  double exactLogZ;
  
  public double getExactLogZ()
  {
    return exactLogZ;
  }
  
  List<Integer> generateData(Random random)
  {
    List<Integer> prior = samplePrior(random);
    List<Integer> result = new ArrayList<>();
    for (int i = 0; i < len; i++)
    {
      int latent = prior.get(i);
      result.add(random.nextCategorical(emissionPrs[latent]));
    }
    return result;
  }
  
  List<Integer> samplePrior(Random random)
  {
    List<Integer> prior = new ArrayList<>();
    int latent = random.nextCategorical(initialPrs);
    prior.add(latent);
    for (int i = 1; i < len; i++)
    {
      latent = random.nextCategorical(transitionPrs[latent]);
      prior.add(latent);
    }
    return prior;
  }
  
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
  
  class HMMKernels implements AnnealingKernels<AnnealedLikelihoodParticle<List<Integer>>>
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
      if (random.nextBernoulli(acceptPr))
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
  
  @Implementations({SimpleThreeStates.class, SimpleTwoStates.class, LowPrValley.class})
  static interface Parameters
  {
    double [][] transitionPrs();
    double [][] emissionPrs();
    double [] initialPrs();
  } 
  
  public static class LowPrValley implements Parameters
  {
    @Arg     @DefaultValue("0.01")
    public double epsilon = 0.01;

    @Override
    public double[][] transitionPrs()
    {
      double [][] result = new double[5][5];
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
          if (Math.abs(i - j) < 2)
            result[i][j] = 1.0;
      for (int r = 0; r < 5; r++)
        Multinomial.normalize(result[r]);
      return result;
    }

    @Override
    public double[][] emissionPrs()
    {
      double [][] result = new double[5][];
      if (epsilon >= 1.0)
        throw new IllegalArgumentException();
      result[0] = result[3] = new double[]{0.2, 0.8};
      result[1] = new double[]{0.1, 0.9};
      result[4] = new double[]{0.3, 0.7};
      result[3] = new double[]{epsilon, 1.0 - epsilon};
      return result;
    }

    @Override
    public double[] initialPrs()
    {
      double [] result = new double[5];
      for (int i = 0; i < 5; i++)
        result[i] = 1.0/5.0;
      return result;
    }
  }
  
  public static class SimpleThreeStates implements Parameters
  {
    @DesignatedConstructor
    public SimpleThreeStates() {}
    
    @Override
    public double[][] transitionPrs()
    {
      return new double[][]{{0.8,0.15,0.05},{0.02,0.93,0.05},{0.15,0.15,0.7}};
    }
    @Override
    public double[][] emissionPrs()
    {
      return transitionPrs();
    }
    @Override
    public double[] initialPrs()
    {
      return new double[]{0.25, 0.25, 0.5};
    }
  }
  
  public static class SimpleTwoStates implements Parameters
  {
    @DesignatedConstructor
    public SimpleTwoStates() {}
    
    @Override
    public double[][] transitionPrs()
    {
      return new double[][]{{0.1, 0.9},{0.6, 0.4}};
    }
    @Override
    public double[][] emissionPrs()
    {
      return transitionPrs();
    }
    @Override
    public double[] initialPrs()
    {
      return new double[]{0.2, 0.8};
    }
  }

}
