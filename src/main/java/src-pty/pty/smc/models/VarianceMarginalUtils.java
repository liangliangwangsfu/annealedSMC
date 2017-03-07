package pty.smc.models;
import java.io.*;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;

import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MathException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math.util.ContinuedFraction;



import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;

import nuts.math.Fct;
import nuts.math.MCIntegrator;
import nuts.math.Sampler;
import nuts.math.Sampling;
import nuts.math.TrapezoidLogSpaceIntegrator;
import nuts.maxent.SloppyMath;
import nuts.util.MathUtils;
import nuts.util.MathUtils.FPlusDelta;

public final class VarianceMarginalUtils
{
  /**
   * Functionalities: approx normalization, sampling and normalized-sampling-density of
   * truncated Generalized Inverse Gaussian distributions
   * 
   * Technique: use closed forms to compute the untruncated normalization,
   * mean and variance (in terms of modified bessel functions of the 
   * second kind, regularized incomplete gamma functions, and std. gamma 
   * and elementary functions), then use moment matching  to
   * fit a gamma distribution (cdf of gamma is a reg. incomp. gamma again)
   * 
   * @author bouchard
   *
   */
//  public static class TruncatedGIG
//  {
//    private final double a, b, p, trunc;
//    
//    
//  }
  
  
//  public static double integrate(
//      final double lambda, 
//      final double oneHalfSumCenteredSquares, 
//      final int n,
//      final double deltaChildHeightPlusVarianceMessages)
//  {
//    final double p = 1.0 - ((double) n) /2.0;
//    final double a = lambda / 2;
//    final double b = oneHalfSumCenteredSquares;
//    final double mean = mean(a,b,p);
//  }
//
//  private static double mean(final double a, final double b, final double p)
//  {
//    final double sqrta = Math.sqrt(a),
//                 sqrtb = Math.sqrt(b);
//    return (sqrtb * K(Bessel.
//  }
  
//  public static double sqrtpiOver2 = Math.sqrt(Math.PI/2);
//  public static double modifiedBesselSecondKind(final double p, final double z)
//  {
//    double pMinusOneHalf = Math.abs(p) - 0.5;
//    if (MathUtils.nint(pMinusOneHalf) - pMinusOneHalf > 0)
//      throw new RuntimeException();
//    final int n = (int) pMinusOneHalf;
//    final double base = sqrtpiOver2 * Math.exp(-z) / Math.sqrt(z);
//    double sum = 0.0;
//    for (int j = 0; j <= n; j++)
//      sum += Math.exp( MathUtils.logFactorial(j + n) 
//          - j * Math.log(2*z) - MathUtils.logFactorial(j) - MathUtils.logFactorial(n-j) );
////    System.out.println(base*sum);
//    System.out.println("std:" + base*sum);
//    System.out.println("logspace:" + Math.exp(logModifiedBesselSecondKind(p,z)));
//    return base * sum;
//  }
  public static double logModifiedBesselSecondKind(final double p, final double z)
  {
    double pMinusOneHalf = Math.abs(p) - 0.5;
    if (MathUtils.nint(pMinusOneHalf) - pMinusOneHalf > 0)
      throw new RuntimeException();
    final int n = (int) pMinusOneHalf;
    final double base = 0.5 * Math.log(Math.PI/2) - z - 0.5*Math.log(z);
    double sum = Double.NEGATIVE_INFINITY;
    for (int j = 0; j <= n; j++)
      sum = SloppyMath.logAdd(sum, ( MathUtils.logFactorial(j + n) 
          - j * Math.log(2*z) - MathUtils.logFactorial(j) - MathUtils.logFactorial(n-j) ));
//    System.out.println(base*sum);
    return base + sum;
  }
//  public static double _modifiedBesselSecondKindOrderNPlusHalf(final int n, final double z)
//  {
//    final double base = sqrtpiOver2 * Math.exp(-z) / Math.sqrt(z);
//    if (n == 0) return base;
//    double current = Double.NaN;
//    double prev = base;
//    double prevPrev = base;
//    for (int i = 1; i <= n; i++)
//    {
//      final double currentNu = i + 0.5;
//      current = prevPrev + (2 * currentNu / z) * prev;
//      prevPrev = prev;
//      prev = current;
//    }
//    return current;
//  }
  
  // test bessel and normalization of GIG 
  public static class GIGlogDensity implements UnivariateRealFunction
  {
    private final double p, a, b;
    public GIGlogDensity(double p, double a, double b)
    {
      this.p = p;
      this.a = a;
      this.b = b;
    }
    public double value(double x) throws FunctionEvaluationException
    {
      if (x <= 0) return Double.NEGATIVE_INFINITY;
      return (p-1)*Math.log(x) + (-(a*x+b/x)/2.0);
//      return Math.pow(x, p-1) * Math.exp(-(a*x+b/x)/2.0);
    }
  }
  
  public static double analyticGIGLogNormalization(final double p, final double a, final double b)
  {
    return Math.log(2) + logModifiedBesselSecondKind(p, Math.sqrt(a*b))
             - (p/2) * Math.log(a/b) ;
  }
  
  public static Pair<Double,Double> meanAndVarOfGIG(double p, double a, double b)
  {
    final double sqrtab = Math.sqrt(a*b);
    final double logKp   = logModifiedBesselSecondKind(p  ,sqrtab),
                 logKpp1 = logModifiedBesselSecondKind(p+1,sqrtab),
                 logKpp2 = logModifiedBesselSecondKind(p+2,sqrtab);
    final double mean = Math.sqrt(b) / Math.sqrt(a) * Math.exp(logKpp1 - logKp);
//    final double mean = Math.sqrt(b) * Kpp1 / Math.sqrt(a) / Kp;
    final double vart2 = Math.exp(logKpp1 - logKp);
//    final double vart2 = Kpp1 / Kp;
    final double var = b/a*(Math.exp(logKpp2-logKp) - vart2*vart2);
//    final double var  = b/a*(Kpp2/Kp - vart2*vart2);
    return Pair.makePair(mean,var);
  }
  
  public static final double TOL = 1e-10;
  public static final int MAX_ITERS = 100000;
//  public static double logRegularizedGammaP(final double a, final double x)
//  {
//    double sum = Double.NEGATIVE_INFINITY;
//    double prev = Double.POSITIVE_INFINITY;
//    for (int n = 0; n < MAX_ITERS; n++) 
//    {
//      sum = SloppyMath.logAdd(sum, n * Math.log(x) - Gamma.logGamma(a+n+1));
//      System.out.println("sum="+sum);
////      System.out.println("prev="+prev);
////      System.out.println("err=" + Math.abs((prev - sum)/sum));
////      System.out.println("weird"+ (Double.POSITIVE_INFINITY < TOL));
//      if (Math.abs((prev - sum)/sum) <  TOL)
//        return sum + a * Math.log(x) - x;
//      prev = sum;
//    }
//    throw new RuntimeException();
//  }
//  public static double logRegularizedGammaQ(final double s, final double z)
//  {
//    final double factor = 0;//s*Math.log(z) - z;
//    double 
//      App = 0,  // Double.NEGATIVE_INFINITY;
//      Ap  = Math.pow(z,s) * Math.exp(-z), ///s*Math.log(z) - z,
//      Bpp = 1, //0,
//      Bp  = (1-s+z); //Math.log(2 - s + z);
//    double previousLogRatio = factor + Math.log(Ap) - Math.log(Bp);
//    for (int i = 2; i < MAX_ITERS; i++) 
//    {
//      final double 
//        a = (i-1)*(s-i), //Math.log(i-1) + Math.log(s-i),
//        b = (2*i-1-s+z); //Math.log(2*i-s+z);
//      double 
//        A = b*Ap + a*App, ///SloppyMath.logAdd(b+Ap, a+App),
//        B = b*Bp + a*Bpp; ///SloppyMath.logAdd(b+Bp, a+Bpp);
//      final double
//        currentLogRatio = factor + Math.log(A) - Math.log(B);
////      System.out.println(A-B);
//      if (Math.abs((currentLogRatio-previousLogRatio)/currentLogRatio) < TOL)
//        return currentLogRatio - Gamma.logGamma(s);
//      App = Ap; Ap = A;
//      Bpp = Bp; Bp = B;
//      previousLogRatio = currentLogRatio;
//    }
//    throw new RuntimeException();
//  }
//  
//  public static double logRegularizedGammaP(final double s, final double z)
//  {
//    double sum = Double.NEGATIVE_INFINITY;
//    double prev = Double.POSITIVE_INFINITY;
//    for (int k = 0; n < 2* MAX_ITERS; n+=2) 
//    {
//      sum = SloppyMath....
//      
//    }
//  }
  
  
  public static double truncatedGIGLogNormalizationApprox(double p, double a, double b, double lowerIntegrand)
  {
    System.out.println("tgiglna("+p+","+a+","+b+","+lowerIntegrand+")");
    if (lowerIntegrand == 0.0)
      return analyticGIGLogNormalization(p,a,b);
    // find mean and var analytically
    final Pair<Double,Double> meanVar = meanAndVarOfGIG(p,a,b);
    final double theta = meanVar.getSecond() / meanVar.getFirst();
    final double kParam = meanVar.getFirst() / theta;
    try
    {
//      System.out.println("rgp(" + kParam + "," + (lowerIntegrand/theta) + ")=" + 
//          logRegularizedGammaQ(kParam, lowerIntegrand/theta));
//      System.out.println(analyticGIGLogNormalization(p,a,b));
      return analyticGIGLogNormalization(p,a,b) +  
        logRegularizedGammaQ(kParam, lowerIntegrand/theta);
    } catch (MathException e) { throw new RuntimeException(e); }
//    return Math.exp(Math.log(analyticGIGNormalization(p,a,b)) 
//        + Math.log(cern.jet.stat.Gamma.incompleteGamma(kParam, lowerIntegrand/theta))
//        - Gamma.logGamma(kParam));
//    return Math.exp(kParam * Math.log(theta) + Gamma.logGamma(kParam));
  }
  
  public static void testGIGNorm(double p, double a, double b, double trunc) throws MaxIterationsExceededException, FunctionEvaluationException, IllegalArgumentException
  {
//    final double analytic = analyticGIGNormalization(p,a,b);
//    System.out.println("Analytic:" + analytic);
    System.out.println("Gamma approx log norm  :" + truncatedGIGLogNormalizationApprox(p,a,b,trunc));
    System.out.println("Numeric approx log norm:" + new TrapezoidLogSpaceIntegrator(new GIGlogDensity(p,a,b)).integrate(trunc, 100));
  }
  
  public static void main(String [] args) throws Exception
  {
    
    {
      double s = 5;
      double x = 10;
      for (int i = 0; i < 1000; i ++)
      {
        System.out.println("x=" + x +",cur:"+ logRegularizedGammaQ(s,x));
        x *= 0.5;
      }
      
    }
    
    if (true) return;
//    -5754.5,5.0,0.07061787167053435,0.38966364702146933
    {
      double 
        p = -5754.5,
        a = 5.0,
        b = 0.07061787167053435,
        trunc = 0.38966364702146933;
      System.out.println("Should converge to:" + truncatedGIGLogNormalizationApprox(p,a,b,0.0));
      for (int i = 0; i < 100; i++)
      {
        System.out.println("current:" + truncatedGIGLogNormalizationApprox(p,a,b,trunc));
        trunc *= 0.5;
      }
      
    }
    if (true) return;
    
    //rgp(5752.499620169006,5.2179022439931566E8)
    double logResult = logRegularizedGammaQ(5752.499620169006, 5.2179022439931566E8);
    System.out.println(logResult);
    
    
    final double a = 5.1, x = 1;
    
    
//    for (int cx = 10; cx < 10000; cx*=2)
//      System.out.println("" + cx + "\t" + Gamma.regularizedGammaP(a, cx));
//    
    System.out.println("Std:" + (Math.log(Gamma.regularizedGammaQ(a, x))));
//    System.out.println("largexapprox:" + (Gamma.logGamma(x)-Gamma.logGamma(a)));
    System.out.println("Log:" + logRegularizedGammaQ(a,x));
//    System.out.println("Log2:" + (1-Math.exp(logRegularizedGammaQ(a,x))));
    //testGIGNorm(-500.5, 3, 300, .05);
//    LogInfo.track("..",false);
//    System.out.println("Numeric:" + (new TrapezoidLogSpaceIntegrator(new UnivariateRealFunction() {
//      public double value(double t) throws FunctionEvaluationException
//      {
//        return -t + (a-1)*Math.log(t);
//      }
//    }).integrate(0,x)-Gamma.logGamma(a)));
//    LogInfo.end_track();
  }
  
  ////////////
  
//  /**
//   * Returns the regularized gamma function P(a, x).
//   * 
//   * @param a the a parameter.
//   * @param x the value.
//   * @return the regularized gamma function P(a, x)
//   * @throws MathException if the algorithm fails to converge.
//   */
//  public static double logRegularizedGammaP(double a, double x)
//      throws MathException
//  {
//      return logRegularizedGammaP(a, x, TOL, Integer.MAX_VALUE);
//  }
//      
//      
//  /**
//   * Returns the regularized gamma function P(a, x).
//   * 
//   * The implementation of this method is based on:
//   * <ul>
//   * <li>
//   * <a href="http://mathworld.wolfram.com/RegularizedGammaFunction.html">
//   * Regularized Gamma Function</a>, equation (1).</li>
//   * <li>
//   * <a href="http://mathworld.wolfram.com/IncompleteGammaFunction.html">
//   * Incomplete Gamma Function</a>, equation (4).</li>
//   * <li>
//   * <a href="http://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html">
//   * Confluent Hypergeometric Function of the First Kind</a>, equation (1).
//   * </li>
//   * </ul>
//   * 
//   * @param a the a parameter.
//   * @param x the value.
//   * @param epsilon When the absolute value of the nth item in the
//   *                series is less than epsilon the approximation ceases
//   *                to calculate further elements in the series.
//   * @param maxIterations Maximum number of "iterations" to complete. 
//   * @return the regularized gamma function P(a, x)
//   * @throws MathException if the algorithm fails to converge.
//   */
//  public static double logRegularizedGammaP(double a, 
//                                         double x, 
//                                         double epsilon, 
//                                         int maxIterations) 
//      throws MathException
//  {
//      double ret;
//
//      if (Double.isNaN(a) || Double.isNaN(x) || (a <= 0.0) || (x < 0.0)) {
//          ret = Double.NaN;
//      } else if (x == 0.0) {
//          ret = 0.0;
//      } else if (a >= 1.0 && x > a) {
//          // use regularizedGammaQ because it should converge faster in this
//          // case.
//          ret = 1.0 - regularizedGammaQ(a, x, epsilon, maxIterations);
//      } else {
//          // calculate series
//          double n = 0.0; // current element index
//          double an = 1.0 / a; // n-th element in the series
//          double sum = an; // partial sum
//          while (Math.abs(an) > epsilon && n < maxIterations) {
//              // compute next element in the series
//              n = n + 1.0;
//              an = an * (x / (a + n));
//
//              // update partial sum
//              sum = sum + an;
//          }
//          if (n >= maxIterations) {
//              throw new MaxIterationsExceededException(maxIterations);
//          } else {
//              ret = Math.exp(-x + (a * Math.log(x)) - logGamma(a)) * sum;
//          }
//      }
//
//      return ret;
//  }
  
  /**
   * Returns the regularized gamma function Q(a, x) = 1 - P(a, x).
   * 
   * @param a the a parameter.
   * @param x the value.
   * @return the regularized gamma function Q(a, x)
   * @throws MathException if the algorithm fails to converge.
   */
  public static double logRegularizedGammaQ(double a, double x)
      throws MathException
  {
      return logRegularizedGammaQ(a, x, TOL, Integer.MAX_VALUE);
  }
  
  /**
   * Returns the regularized gamma function Q(a, x) = 1 - P(a, x).
   * 
   * The implementation of this method is based on:
   * <ul>
   * <li>
   * <a href="http://mathworld.wolfram.com/RegularizedGammaFunction.html">
   * Regularized Gamma Function</a>, equation (1).</li>
   * <li>
   * <a href="    http://functions.wolfram.com/GammaBetaErf/GammaRegularized/10/0003/">
   * Regularized incomplete gamma function: Continued fraction representations  (formula 06.08.10.0003)</a></li>
   * </ul>
   * 
   * @param a the a parameter.
   * @param x the value.
   * @param epsilon When the absolute value of the nth item in the
   *                series is less than epsilon the approximation ceases
   *                to calculate further elements in the series.
   * @param maxIterations Maximum number of "iterations" to complete. 
   * @return the regularized gamma function P(a, x)
   * @throws MathException if the algorithm fails to converge.
   */
  public static double logRegularizedGammaQ(final double a, 
                                         double x, 
                                         double epsilon, 
                                         int maxIterations) 
      throws MathException
  {
      double ret;

      if (Double.isNaN(a) || Double.isNaN(x) || (a <= 0.0) || (x < 0.0)) {
          ret = Double.NaN;
      } else if (x == 0.0) {
          ret = 0.0; //1.0;
//      } else if (x < a || a < 1.0) {
//          // use regularizedGammaP because it should converge faster in this
//          // case.
//          ret = 1.0 - regularizedGammaP(a, x, epsilon, maxIterations);
      } else {
          // create continued fraction
          ContinuedFraction cf = new ContinuedFraction() {

              private static final long serialVersionUID = 5378525034886164398L;

              protected double getA(int n, double x) {
                  return ((2.0 * n) + 1.0) - a + x;
              }

              protected double getB(int n, double x) {
                  return n * (a - n);
              }
          };
          
          //ret = 1.0 / cf.evaluate(x, epsilon, maxIterations);
          ret = -Math.log(cf.evaluate(x, epsilon, maxIterations));
//          ret = Math.exp(-x + (a * Math.log(x)) - logGamma(a)) * ret;
          ret = (-x + (a * Math.log(x)) - Gamma.logGamma(a)) + ret;
      }

      return ret;
  }

}
