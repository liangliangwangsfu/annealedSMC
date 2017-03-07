package pty.smc;
import pty.eval.SymmetricDiff;
import pty.learn.DiscreteBP;
import goblin.Taxon;
import java.io.*;
import java.util.*;

import conifer.particle.PhyloParticle;

import monaco.prop.Proposal;
import nuts.math.GMFct;
import nuts.math.Sampling;
import nuts.math.TreeSumProd;
import nuts.maxent.SloppyMath;
import nuts.util.Counter;
import nuts.util.MathUtils;

import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;

import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.models.CTMC;

/**
 * See Yee Whye Teh et al. Bayesian Agglomerative Clustering with Coalescents
 * @author bouchard
 */
public  class PriorPriorKernel2 implements LazyParticleKernel<PartialCoalescentState[]>
{
  @Option public static boolean printBranchLengthMagnitudes = false;
//  @Option public static double [] r=new double[]{0.02121, 0.15549, 0.49708, 1.10712, 3.24910}; 
  public double [] r=null; //new double[]{0.02121, 0.15549, 0.49708, 1.10712, 3.24910};
//  @Option public  double [] r=new double[]{1};
  private final PartialCoalescentState [] initial;
  public PriorPriorKernel2(PartialCoalescentState[] initial) { this.initial = initial; }
  public PartialCoalescentState[] getInitial() { return initial; }
    
  @Override
  public int nIterationsLeft(PartialCoalescentState[] partialState) {
	    return partialState[0].nIterationsLeft();
	    }

  
  public Object _next(
      Random rand,
      PartialCoalescentState[] current, boolean isPeek)
  {	
	  
	  PartialCoalescentState[] result = new PartialCoalescentState[r.length];
	  double [] wVec=new double[r.length], oldwVec=new double[r.length];
	  double w = 0, oldw=0;
	    // 1- sample the exponential waiting time
//	    final double delta = Sampling.sampleExponential(rand, (initial.isClock() ? 1.0 : 0.5)/nChoose2(current.nRoots()));        
	  final double delta = Sampling.sampleExponential(rand, (initial[0].isClock() ? 1.0 : 0.5)/nChoose2(current[0].nRoots()));
//		  final double delta = Sampling.sampleExponential(rand, (initial[0].isClock() ? 0.1 : 0.05)/nChoose2(current[0].nRoots()));
	 //   System.out.println("rand: "+rand.toString()+"	rate: "+(initial.isRooted() ? 1.0 : 0.5)/nChoose2(current.nRoots()) +"	delta is: "+delta);
	    
	    // 2- sample a random pair (without replacement)
	    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, current[0].nRoots(), 2);
	    
	    final int i0 = sampledIndices.get(0),
	              i1 = sampledIndices.get(1);
	    
	    double leftIncrement = 0.0,
	           rightIncrement= 0.0;
	    if (!initial[0].isClock())
	    {
	      //double incr = Sampling.sampleExponential(rand, 0.1);
	    	double incr = Sampling.sampleExponential(rand, 1);
	      if (rand.nextBoolean())
	        leftIncrement += incr;
	      else
	        rightIncrement+= incr;
	    }
	    
	    
	    if (printBranchLengthMagnitudes)
	    {
	      LogInfo.logsForce("Merged things:" + current[0].clade(i0)+ "," + current[0].clade(i1));
	      
	      // show likelihood for different times
	      for (int i = -20; i < 20; i += 2)
	      {
	        PartialCoalescentState result0 = current[0].coalesce(
	            i0, i1, 
	            delta * Math.pow(2.0, i), leftIncrement, rightIncrement);
	        LogInfo.logsForce("\t"+ i + "\t"+ delta * Math.pow(2.0, i) + "\t" + result0.logLikelihoodRatio());
	      }
	    }

	  
	  for(int j=0;j<r.length;j++){

//    if (returnProposalLogRatio)
//    {
//      PhyloParticle result = null;
//      if (!isPeek)
//        result = current.coalesce(
//            i0, i1, 
//            delta, leftIncrement, rightIncrement);
//      return Pair.makePair(result, second)
//    }
            
    if (isPeek)
    {
    	double temp1=current[j].peekLogLikelihoodRatio(i0, i1, delta*r[j], leftIncrement*r[j], rightIncrement*r[j]);
    	double temp2=current[j].logLikelihood();    	         	
    	wVec[j]=current[j].peekLogLikelihoodRatio(i0, i1, delta*r[j], leftIncrement*r[j], rightIncrement*r[j])+current[j].logLikelihood(); 
    	oldwVec[j]=current[j].logLikelihood(); 
//      w = SloppyMath.logAdd(w,current[j].peekLogLikelihoodRatio(i0, i1, delta*r[j], leftIncrement*r[j], rightIncrement*r[j])+current[j].logLikelihood());
//      oldw = SloppyMath.logAdd(oldw, current[j].logLikelihood());
    }
    else
      result[j] = current[j].coalesce(i0, i1, delta*r[j], leftIncrement*r[j], rightIncrement*r[j]);    
	  }
	    
	    // 3- the weight update is simply equal to the ratio of the new likelihood score to the old one
     w=SloppyMath.logAdd(wVec);
     oldw=SloppyMath.logAdd(oldwVec);
	  
	    if (isPeek)
	    {
//	    	System.out.println("w - oldw = "+(w - oldw));
	      return w - oldw;
	    }
	    else{
	    	double [] LikeR=new double[r.length], oldLikeR=new double[r.length];
	    	
	    	for(int k=0;k<r.length;k++){ 
	    		LikeR[k]=result[k].logLikelihood(); 
	    		oldLikeR[k]=current[k].logLikelihood(); 
//	    	     LikeR=SloppyMath.logAdd(LikeR,);
//	    		oldLikeR=SloppyMath.logAdd(oldLikeR,);	    		
//        System.out.println("k "+k +" r[k] "+r[k]+" LikeR[k] "+LikeR[k] +" oldLikeR[k] "+oldLikeR[k]);
	    	}	    	
//	    	System.out.println("SloppyMath.logAdd(LikeR)-SloppyMath.logAdd(oldLikeR) = "+(SloppyMath.logAdd(LikeR)-SloppyMath.logAdd(oldLikeR)));
	      return Pair.makePair(result, SloppyMath.logAdd(LikeR)-SloppyMath.logAdd(oldLikeR));
	    }
	  
	  
  }
  public static double nChoose2(double n) { return n*(n-1)/2; }
  @Override
  public Pair<PartialCoalescentState[], Double> next(Random rand,
      PartialCoalescentState[] current)
  {
    return (Pair) _next(rand, current, false);
  }
  @Override
  public double peekNext(Random rand, PartialCoalescentState[] current)
  {
    return (Double) _next(rand, current, true);
  }

}
