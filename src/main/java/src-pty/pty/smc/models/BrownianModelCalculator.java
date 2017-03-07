package pty.smc.models;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import nuts.util.Counter;

import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.special.Gamma;

import fig.basic.Option;
import fig.basic.Pair;
import fig.prob.Gaussian;
import fig.prob.SampleUtils;
import goblin.Taxon;

import pty.smc.ParticleFilter;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.models.VarianceMarginalUtils.GIGlogDensity;
import pty.smc.test.TestBrownianModel;

public class BrownianModelCalculator implements LikelihoodModelCalculator {
  
  @Option public static boolean useVarianceTransform = true;
  public final boolean resampleRoot;
  
  // sites X {mean,variance}
  public final double[] message;
  public final double messageVariance;
  private final double loglikelihood;
  public double lognorm;
  public final BrownianModel bm;
  private final int nsites;
//  public final double varianceScale;
  
  
  private BrownianModelCalculator (double[] message, double messageVariance, BrownianModel bm, double loglikelihood, boolean resampleRoot) {
    this.message = message;
    this.messageVariance = messageVariance;
    this.loglikelihood = loglikelihood;
    this.bm = bm;
    this.nsites = bm.nsites;
    this.lognorm = 0;
    this.resampleRoot = resampleRoot;
//    this.varianceScale = bm.varianceScale;
  }
  
  
  private BrownianModelCalculator (double[] message, double messageVariance, BrownianModel bm, double loglikelihood, double lognorm, boolean resampleRoot) {
	    this.message = message;
	    this.messageVariance = messageVariance;
	    this.loglikelihood = loglikelihood;
	    this.bm = bm;
	    this.nsites = bm.nsites;
	    this.lognorm = lognorm;
	    this.resampleRoot = resampleRoot;
//	    this.varianceScale = bm.varianceScale;
	  }
  
  public static BrownianModelCalculator observation (double[] observed, BrownianModel bm, boolean resampleRoot) {
    return observation(observed, bm, useVarianceTransform, resampleRoot);
  }
  private static BrownianModelCalculator observation (double[] observed, BrownianModel bm, boolean useVarianceTransform, boolean resampleRoot) {
    if (observed.length !=bm.nsites )
      throw new RuntimeException();
    double[] state = getNewMessage (bm.nsites);
    for (int i = 0 ; i < bm.nsites; i++) {
     if (useVarianceTransform) // Delta method !!
        state[i] = Math.asin(Math.sqrt(observed[i]));
      else // assumes perfect observations
        state[i] = observed[i];
    }
    return new BrownianModelCalculator (state,0.0, bm, 0, resampleRoot);
  }

  private final Object calculate(
      final LikelihoodModelCalculator node1, final LikelihoodModelCalculator node2, 
      final double v1, final double v2,
      final boolean peek) 
  {
//    System.out.println("calculate("+v1+","+v2+")");
    final BrownianModelCalculator 
      l1 = (BrownianModelCalculator) node1,
      l2 = (BrownianModelCalculator) node2;
    final double[] message  = (peek ? null : getNewMessage(nsites));
    double logl = 0;
    final double var1 = l1.messageVariance; // l1.message[i][1];
    final double var2 = l2.messageVariance; //l2.message[i][1];
    final double var = 1/(var1+v1) + 1/(var2+v2);
    final double newMessageVariance = 1/var;
    for (int i = 0 ; i < nsites; i++) 
    {
      final double mean1 = l1.message[i];
      final double mean2 = l2.message[i];

      
      if (!peek)
        message[i] = ( (mean1)/(var1+v1) + (mean2)/(var2+v2) ) / var;
      
      final double cur  = logNormalDensity (mean1-mean2,0,bm.varianceScale*(v1+var1+v2+var2));
      
      logl += cur;
    }
//    logl += l1.loglikelihood + l2.loglikelihood;
    double lognorm = logl;
    logl += l1.loglikelihood + l2.loglikelihood; 
    //(l1.lognorm + l2.lognorm);
//    System.out.println("lognorm="+lognorm+",log="+logl);
    
    if (peek) return logl;
    else      return (resampleRoot ? resampleRoot(message, newMessageVariance, bm) : new BrownianModelCalculator(message, newMessageVariance,  bm, logl, lognorm, resampleRoot)); 
  }
  private static Random rand = new Random(1);
  private static BrownianModelCalculator resampleRoot(double[] message,
      double newMessageVariance, BrownianModel bm)
  {
    for (int i = 0; i < message.length; i++)
    {
      message[i] = Gaussian.sample(rand, message[i], newMessageVariance);
    }
    return observation(message, bm, false, true);
  }


  public LikelihoodModelCalculator combine(LikelihoodModelCalculator node1,
      LikelihoodModelCalculator node2, double delta1, double delta2, boolean doNotBuildCache)
  {
    return (LikelihoodModelCalculator) calculate(node1, node2, delta1, delta2, false);
  }
  
  public final Pair<Double, Double> sampleBranchLength (final CoalescentNode node1, final CoalescentNode node2, double topHeight,  
      int nroots, double particleweight, Random rand) {
    final BrownianModelCalculator l1 = (BrownianModelCalculator) (node1.likelihoodModelCache),
      l2 = (BrownianModelCalculator) (node2.likelihoodModelCache);
    final double   var1 = l1.messageVariance; //l1.variancemessage;
    final double var2 = l2.messageVariance; // l2.variancemessage;
      
    double c = var1 + var2 + 2 * topHeight - node1.height - node2.height;  
    double b = 0;
    double a =  0.5*nroots * (nroots - 1);
    int m = nsites;
    double p = 0.5*(m-2);
    
    for (int i = 0 ; i < nsites; i++) {
      final double mean1 = l1.message[i]; //l1.meanmessage[i];
      final double mean2 = l2.message[i]; //l2.meanmessage[i];
        
   // TODO: INVESTIGATE IF WE CAN ACTUALLY HAVE TRACTABLE SITE-SPECIFIC VARIANCES...
        b += Math.pow(mean1 - mean2, 2)/bm.varianceScale; //variances[i];
    }
    
    
  /*  
    System.out.println ("Sampling from GIG with " 
          + "\tnroots = " + nroots + "\ta = " + a
        + "\tb = " + b
        + "\tc = " + c
        + "\tm = " + m);
*/
    
    try {
    	double mode  = (m-4)/b;
    	double result = 0;
    	double correct= 0;
    	
    	if ( mode < 1/c ) {
    	
    		while (true) {
    			int maxiter = (int)(0.5*m - 2);
      
    			double sum = 0 ;
    			for (int i = 0 ; i < maxiter; i ++)
    				sum -= Math.log(rand.nextDouble());
      
    			GammaDistributionImpl gamma = new GammaDistributionImpl(0.5*m-1 - maxiter, 1);      
    			sum +=  gamma.inverseCumulativeProbability(rand.nextDouble());
      
    			result = sum * (2/b);
    			if (result <= 1/c)
    				break;
    		}
    	    // We are sampling the inverse of the branch length
    	    // because we did a transform on the GIG
//    	    System.out.println( "Sampled  = " + result);  
    	    result = 0.5 * (1/result - c);
    	    correct = -0.5 * m  * (Math.log(2* Math.PI) + Math.log(bm.varianceScale));
    	} else {
//    		System.out.println ("Using uniform proposal");
    		result = rand.nextDouble() * (1/c);

    		GIGlogDensity gig = new GIGlogDensity(p - 1, b, a);
    		
    		BrownianModelCalculator bmc1 = (BrownianModelCalculator) node1.likelihoodModelCache;
    	    BrownianModelCalculator bmc2 = (BrownianModelCalculator) node2.likelihoodModelCache;
    	    correct = gig.value(result) + 0.5*a*c - particleweight -  Math.log(2*c);
   	     	correct += (bmc1.lognorm + bmc2.lognorm - bmc1.logLikelihood() - bmc2.logLikelihood());
    		correct -= (0.5  * m * (Math.log(2*Math.PI) + Math.log(bm.varianceScale)));
    		
//    	    System.out.println( "Sampled  = " + result);  
    	    result = 0.5 * (1/result - c);

    	}
      
 //     System.out.println( "Sampled bl = " + result);
      return Pair.makePair(result, correct);
      
    } catch (Exception e ) {
 //     System.out.println ("Using normal approximation");
      double correct = 0 ;
      double result  = rand.nextGaussian();
      result = (result + (m-2)/b) *(Math.sqrt(2*(m-2))/b);
      
      result = 0.5 * (1/result - c);
//      System.out.println( "Sampled = " + result);
      return Pair.makePair (result, correct);
    }
  }
  
  
  public final double evaluatePair (final CoalescentNode node1, final CoalescentNode node2, 
      int nroots, double topHeight) {
    final BrownianModelCalculator l1 = (BrownianModelCalculator) (node1.likelihoodModelCache),
      l2 = (BrownianModelCalculator) (node2.likelihoodModelCache);
    final double   var1 = l1.messageVariance; //l1.variancemessage;
      final double var2 = l2.messageVariance; // l2.variancemessage;
      
//      System.out.println ("Evaluating pair " + node1.nodeIdentifier + "," + node2.nodeIdentifier);
//      System.out.println ("var1 = " + var1 + "\tvar2 = " + var2 
//                  + "\th1 = " + node1.height + "\t h2 = " + node2.height);
      
    double c = var1 + var2 + 2 * topHeight - node1.height - node2.height;
    double b = 0;
    double a =  0.5 * nroots * (nroots - 1);
    int m = nsites;
    double p = 0.5*(m-2);
    
    for (int i = 0 ; i < nsites; i++) {
        final double mean1 = l1.message[i]; //l1.meanmessage[i];
        final double mean2 = l2.message[i]; //l2.meanmessage[i];
        
        // TODO: INVESTIGATE IF WE CAN ACTUALLY HAVE TRACTABLE SITE-SPECIFIC VARIANCES...
        b += Math.pow(mean1 - mean2, 2)/bm.varianceScale; //variances[i];
    }
    
    
//    System.out.println ("Computing lognorm with " 
//                + "nroots = " + nroots + "\ta = " + a
//                + "\tb = " + b + "\tc = " + c 
//                + "\tm = " + m);
    
    double weight;
    
    if (TestBrownianModel.test) { 
    	double hi = node1.height;
    	double hj = node2.height;
    	double v = (topHeight + 1 - hi)+(topHeight+1-hj);
    	v +=(var1 + var2);
    	weight = -0.5 * Math.log(2 * Math.PI * v);
//    	System.out.println ("nroots = " + nroots);
    	weight -= Math.log(0.5*nroots*(nroots-1));
    		
    	BrownianModelCalculator bmc1 = (BrownianModelCalculator) node1.likelihoodModelCache;
        BrownianModelCalculator bmc2 = (BrownianModelCalculator) node2.likelihoodModelCache;
        weight += (bmc1.lognorm + bmc2.lognorm - bmc1.logLikelihood() - bmc2.logLikelihood());
        
        if (Double.isInfinite(weight) || Double.isNaN(weight)) {
//    		System.out.println ("Illegal weight 2");
    		System.exit (1);
    	}
    	return weight;
    	
    	
    	
    }
    
    // Compute the integral
    if (node1.isLeaf() && node2.isLeaf()) {
      weight = GIG.GIGapproxLognorm(b, a, p);
    } else {
      weight = GIG.GIGapproxLognorm(b, a, p);
      if (Double.isNaN(weight)) {
//    	  System.out.println ("weight1 = " + weight);
    	  System.exit (1);
      }
      
      try {
        GammaDistributionImpl gamma = new GammaDistributionImpl(p, 1/(0.5*b));
        double result = Math.log (gamma.cumulativeProbability (0,1/c));
        if (Double.isInfinite(result) || Double.isNaN (result) ) {
          result = logGammaDensity(1/c, p, 0.5*b) + (1/c);
        }

        weight += result;

        if (Double.isNaN(weight)) {
//      	  System.out.println ("weight2 = " + weight);
//      	  System.out.println ("result1 = " + Math.log (gamma.cumulativeProbability (0,1/c)));
//      	  System.out.println ("result2 = " + logGammaDensity (1/c,p,0.5*b));
      	System.exit (1);
        }
        
      } catch (Exception e ) {
        e.printStackTrace ();
        System.exit (1);
      }
    }

    weight += (a*c)/2;
    BrownianModelCalculator bmc1 = (BrownianModelCalculator) node1.likelihoodModelCache;
    BrownianModelCalculator bmc2 = (BrownianModelCalculator) node2.likelihoodModelCache;
    weight += (bmc1.lognorm + bmc2.lognorm - bmc1.logLikelihood() - bmc2.logLikelihood());
    
    if (Double.isNaN(weight)) {
//  	  System.out.println ("weight3 = " + weight);
  	  System.exit (1);
    }
    
    
//    System.out.println ("weight = " + weight);
    return weight;    
  }
  
  public double logGammaDensity (double x, double a, double b ) {
   double y = (a-1)*Math.log(x) - b*x + a * Math.log(b) - Gamma.logGamma(a);
   return y;
  }

  public double peekCoalescedLogLikelihood(LikelihoodModelCalculator node1,
      LikelihoodModelCalculator node2, double delta1, double delta2)
  {
    return (Double) calculate(node1, node2, delta1, delta2, true);
  }
  
  public static final double logNormalDensity (final double x, final double mean, final double var) {
    return -0.5*(x-mean)*(x-mean)/var -0.5*Math.log(2*Math.PI * var);
  }
  
  public double extendLogLikelihood(double delta) {
    return loglikelihood;
  }

  public double logLikelihood() {
    return loglikelihood;
  }
  
  private static double[] getNewMessage (final int nsites) {
    return new double[nsites];
  }
//    if (masterCacheEnabled && threadSpecificCacheEnabled())
//    {
//      double [][] result = null;
//      int threadIndex = threadMap.get(Thread.currentThread());
//      result = cacheCache[threadIndex];
//      if (result != null)
//        return result;
//      else
//        return cacheCache[threadIndex] = new double[nsites][2];
//    }
//    return new double[nsites][2];
////    for (int i = 0 ; i < nsites; i++)
////      state[i] = new double[2];
////    return state;
//  }
//  
//  private static boolean threadSpecificCacheEnabled()
//  {
//    return cacheCacheEnabled[threadMap.get(Thread.currentThread())];
//  }
//  public static void setTreadSpecificCache(boolean value)
//  {
//    if (cacheCache == null)
//      cacheCache = new double[nThreads][][];
//    if (value == true) 
//      masterCacheEnabled = true;
//    if (cacheCacheEnabled == null)
//      cacheCacheEnabled = new boolean[100];
//    Integer threadIndex = threadMap.get(Thread.currentThread());
//    if (threadIndex == null)
//      synchronized(BrownianModelCalculator.class)
//      {
//        if (threadMap.size() == nThreads)
//          threadMap = new HashMap<Thread,Integer>();
//        threadIndex = nextIndex(threadMap);
//        threadMap.put(Thread.currentThread(), threadIndex);
//      }
//    cacheCacheEnabled[threadIndex] = value;
//  }
//  
//  private static Integer nextIndex(Map<Thread, Integer> threadMap)
//  {
//    int max = -1;
//    for (int value : threadMap.values())
//      if (value > max)
//        max = value;
//    return max+1;
//  }
//
//  private static double[][][] cacheCache = null;
//  private static Map<Thread,Integer> threadMap = new HashMap<Thread,Integer>();
//  private static boolean[] cacheCacheEnabled = null;
//  private static boolean masterCacheEnabled = false;
//  private static Integer nThreads = null;
//  public static void setNThreads(int n) { nThreads = n; }

  public boolean isReversible() { return true; }
}
