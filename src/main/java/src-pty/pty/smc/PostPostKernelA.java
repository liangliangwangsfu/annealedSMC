package pty.smc;
import java.util.Random;

import nuts.math.Sampling;
import nuts.math.TrapezoidLogSpaceIntegrator;
import nuts.util.CoordinatesPacker;

import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;

import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import pty.smc.models.VarianceMarginalUtils;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Pair;
import fig.prob.SampleUtils;

public class PostPostKernelA implements ParticleKernel<PartialCoalescentState>
{
  private static final double MAX_DELTA = 3;
  private static final int MAX_TRIALS = 1000000;
  private final PartialCoalescentState initial;
  public PostPostKernelA(
      PartialCoalescentState initial) { 
    this.initial = initial;
  }
  public Pair<PartialCoalescentState, Double> next (Random rand, PartialCoalescentState state) {
    int nRoots = state.nRoots();
    double [] prs = new double[nRoots*nRoots];
    double lognorm = 0.0;
    CoordinatesPacker cp = new CoordinatesPacker(nRoots);
    for (int i = 0; i < nRoots; i++)
      for (int j = 0; j < nRoots; j++)
      {
        final int idx = cp.coord2int(i,j);
        if (i >= j) prs[idx] = Double.NEGATIVE_INFINITY;
        else        lognorm += prs[idx] = posteriorCoalesceLogPr(i,j,state); 
      }
    NumUtils.expNormalize(prs);
    final int idx = SampleUtils.sampleMultinomial(rand, prs); 
    final int [] pair = cp.int2coord(idx);
    final double delta = Sampling.sampleExponential(rand, (initial.isClock() ? 1.0 : 0.5)/PriorPriorKernel.nChoose2(state.nRoots())); //samplePosteriorDelta(rand, state, pair[0], pair[1]);   
    PartialCoalescentState result = state.coalesce(pair[0],pair[1],delta,0,0); // USE ROOTED ASSUMPTION HERE
    final double weightUpdate = prs[idx] + lognorm - result.logLikelihoodRatio();
    return Pair.makePair(result, weightUpdate);
  }
  private double samplePosteriorDelta(Random rand,
      PartialCoalescentState state, int left, int right)
  {
    if (!state.isBrownianMotion())
      throw new RuntimeException(); // need something else for CTMC
    Pair<UnivariateRealFunction, Double> pair 
      = getProposal(left, right, state);
    final double M = pair.getSecond();
    final double param = 1.0 / PriorPriorKernel.nChoose2(state.nRoots());
    LogInfo.track("Sampling");
    for (int i = 0; i < MAX_TRIALS; i++) 
    {
      LogInfo.logs("Attempt " + i);
      final double x = Sampling.sampleExponential(rand, param);
      final double u = rand.nextDouble() * M;
      try {
        if (u <= pair.getFirst().value(x))
        {
          LogInfo.end_track();
          return x;
        }
      } catch (FunctionEvaluationException e) { throw new RuntimeException(); }
    }
    throw new RuntimeException(); // TODO: fix this and sample from prior or something
  }
  /**
   * @param left
   * @param right
   * @param state
   * @return <proposal, upper bound on the density>
   */
  private static Pair<UnivariateRealFunction,Double> getProposal(int left, int right,
      final PartialCoalescentState state)
  {
    final LikelihoodModelCalculator 
      leftCalculator  = state.getLikelihoodModelCalculator(left),
      rightCalculator = state.getLikelihoodModelCalculator(right);
    final double  
      node1H = state.getHeight(left),
      node2H = state.getHeight(right);
    final double  // compute  baseV1, baseV2; USE ROOTED ASSUMPTION HERE
      baseV1 = Math.max(0,node2H - node1H),
      baseV2 = Math.max(0,node1H - node2H);
    final BrownianModelCalculator 
    l1 = (BrownianModelCalculator) state.getLikelihoodModelCalculator(left),
    l2 = (BrownianModelCalculator) state.getLikelihoodModelCalculator(right);
    final double halfSumSquareDiff = PostPostKernelA.halfSumSquareDiff(l1.message, l2.message);
    final double n = state.getObservations().nSites();
    final double prefix = (n/2) * Math.log(2*Math.PI);
    // create integrand
    UnivariateRealFunction fct = new UnivariateRealFunction() {
      public double value(double deltaT) throws FunctionEvaluationException
      {
        
        // zero creates NaN?
        if (baseV1 + deltaT + baseV2 + deltaT == 0) deltaT = 1e-10;
        final double sigma = 2*deltaT + Math.abs(node1H-node2H) + l1.messageVariance + l2.messageVariance;
//        System.out.println("sigma"+ sigma);
        final double secondTerm = -prefix - (n/2) * Math.log(sigma) - halfSumSquareDiff / sigma;
//        System.out.println("1=" + secondTerm);
//        System.out.println("2=" + leftCalculator.peekCoalescedLogLikelihood(
//                leftCalculator, rightCalculator, baseV1 + deltaT, baseV2 + deltaT, false));
        final double result = logPrior(deltaT, state.nRoots()) + secondTerm;
            
//            leftCalculator.peekCoalescedLogLikelihood(
//                leftCalculator, rightCalculator, baseV1 + deltaT, baseV2 + deltaT, false);
//        if (result == 0.0)
//          throw new RuntimeException("Probably underflow; before exp, value was:" + 
//              (logPrior(deltaT, state.nRoots()) + 
//              leftCalculator.peekCoalescedLogLikelihood(
//                  leftCalculator, rightCalculator, baseV1 + deltaT, baseV2 + deltaT, false)));
        return result;
      }
    };
    return Pair.makePair(fct, leftCalculator.peekCoalescedLogLikelihood(leftCalculator, rightCalculator, baseV1, baseV2));
  }
  private static double posteriorCoalesceLogPr(int left, int right,
      PartialCoalescentState state)
  { 
    if (!state.isReversible() || !state.isClock())
      throw new RuntimeException("Not yet supported");
    
    
    
    double sum = 0;
    // first, add all the log pr of things that do not coalesce 
    // USE REVERSIBILITY HERE
//    for (int i = 0; i < state.nRoots(); i++)
//      if (i!=left && i!=right)
//        sum += state.getLikelihoodModelCalculator(i).logLikelihood();
    // THIS CANCELS OUT b/c when we merge we need to add the likelihood of the two leaves anyways...
    try {
      LogInfo.track("Integrating...",true);
      TrapezoidLogSpaceIntegrator tlsi = new TrapezoidLogSpaceIntegrator(getProposal(left,right,state).getFirst());
//      tlsi.setRelativeAccuracy(0.1);
      sum += tlsi.integrate(0, MAX_DELTA);
      LogInfo.end_track();
    } catch (Exception e) { throw new RuntimeException(e); }
    /////// 
/////// 
    // compare to gamma approximation
    {
      System.out.println("Numerical approx.  :" + sum);
      final BrownianModelCalculator 
        bmc1 = (BrownianModelCalculator) state.getLikelihoodModelCalculator(left),
        bmc2 = (BrownianModelCalculator) state.getLikelihoodModelCalculator(right);
      final double 
        h1 = state.getHeight(left),
        h2 = state.getHeight(right);
      final double lambda = 1.0/PriorPriorKernel.nChoose2(state.nRoots());
      final double k = Math.abs(h1-h2) + bmc1.bm.varianceScale * (bmc1.messageVariance + bmc2.messageVariance);
      final double kPrime = halfSumSquareDiff(bmc1.message,bmc2.message);
      final double n = state.getObservations().nSites();
      final double p = 1.0 - n/2;
      final double a = lambda/2/bmc1.bm.varianceScale;
      final double truncIntegral = VarianceMarginalUtils.truncatedGIGLogNormalizationApprox(p,a,kPrime, k);
      final double estimate = 
        Math.log(lambda) - Math.log(2*bmc1.bm.varianceScale) + k/2/bmc1.bm.varianceScale - (n/2) * Math.log(2*Math.PI) + truncIntegral;
     
      System.out.println("Gamma approximation:" + estimate);
    }
    //
    //////
    
    
    //////
    return sum;
  }
  private static double halfSumSquareDiff(double[] message1, double[] message2)
  {
    if (message1.length != message2.length) throw new RuntimeException();
    double sum = 0;
    for (int i =0 ;i < message1.length; i++)
    {
      if (message1[i]<0 || message1[i]>Math.PI/2) throw new RuntimeException();
      if (message2[i]<0 || message2[i]>Math.PI/2) throw new RuntimeException();
      final double term = message1[i] - message2[i];
      sum += term*term;
    }
    return sum/2;
  }
  public static double logPrior(double x, int nRoots)
  {
    final double param = 1.0/PriorPriorKernel.nChoose2(nRoots);
    return Sampling.exponentialDensity(param,x);
  }
  public PartialCoalescentState getInitial() { return initial; }
  
  public int nIterationsLeft(PartialCoalescentState partialState)
  {
    return partialState.nIterationsLeft();
  }
}
