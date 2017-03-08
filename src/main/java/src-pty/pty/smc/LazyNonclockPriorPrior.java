package pty.smc;

//import ps.LazyPS;
import pty.eval.SymmetricDiff;
import pty.learn.DiscreteBP;
import goblin.Taxon;
import java.io.*;
import java.util.*;

import nuts.math.GMFct;
import nuts.math.Sampling;
import nuts.math.TreeSumProd;
import nuts.util.Arbre;
import nuts.util.Counter;
import nuts.util.MathUtils;

import ev.ex.PathCount;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;

import pty.smc.LazyPCS;
import pty.smc.LazyPS;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleKernel;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.models.CTMC;

/**
 */
public class LazyNonclockPriorPrior implements ParticleKernel<LazyPCS>
{
  @Option public boolean useCountCorrection=true;
  @Option public double nonclockTreeRate=1.0;	
  private final LazyPCS initial;
  public LazyNonclockPriorPrior(PartialCoalescentState initial) { this.initial = new LazyPCS(initial); }
  public LazyPCS getInitial() { return initial; }
  public int nIterationsLeft(LazyPCS partialState)
  {
    return partialState.nIterationsLeft();
  }
  
  public Pair<LazyPCS, Double> next(
      Random rand,
      LazyPCS _current)
  {
    PartialCoalescentState current = _current.getState();
    int forestSize = current.nRoots();
    // 1- sample the two branch lengths
//    double localTreeRate=nonclockTreeRate*forestSize;
//    double localTreeRate=nonclockTreeRate*nChoose2(forestSize);
    double localTreeRate=nonclockTreeRate;
    double leftBranchLength = Sampling.sampleExponential(rand, 1.0/localTreeRate),
    		rightBranchLength=Sampling.sampleExponential(rand, 1.0/localTreeRate);
    
    // 2- sample a random pair (without replacement)
    List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, forestSize, 2);
    
    final int i0 = sampledIndices.get(0),
              i1 = sampledIndices.get(1);
    
//    System.out.println("index ("+i0 +", "+i1+") "+ "lengths: ("+leftBranchLength+", "+rightBranchLength+")");
    double logll = current.peekLogLikelihoodRatio(i0, i1, 0, leftBranchLength,rightBranchLength);
//   System.out.println("Branch lengths: "+leftBranchLength+", "+rightBranchLength+ " loglikelihood ratio: "+logll); 
 
    double logCountRatio=0;
    if(useCountCorrection) logCountRatio=logCountRatio(current, i0, i1);
        
    // 3- the weight update is simply equal to the ratio of the new likelihood score to the old one
    return Pair.makePair(new LazyPCS(current, 0, i0, i1, leftBranchLength,rightBranchLength), logll+logCountRatio);
//    return Pair.makePair(new LazyPS(current, i0, i1, leftBranchLength,rightBranchLength), logll);
  }
  
  
//  public static double logCount(PartialCoalescentState current)
//  {
//	  if(current.nRoots()==2)  
//		  {
//		  Arbre<Taxon> firstArbre=current.getSubtree(0).topology(), secondArbre=current.getSubtree(1).topology();
//		  long first=MathUtils.nint(PathCount.pathLength(firstArbre)),
//		  second=MathUtils.nint(PathCount.pathLength(secondArbre));
////System.out.println(PathCount.logChoose((int)(first+second), (int) first)+" "+PathCount.logPathCount(firstArbre)+" "+PathCount.logPathCount(secondArbre));
//		  return PathCount.logChoose((int)(first+second), (int) first)+PathCount.logPathCount(firstArbre)+PathCount.logPathCount(secondArbre);
//		  }
//	  
//	  long [] countVec=new long[current.nRoots()];
//	  long total=0;
//	  double result=0;
//	  for(int i=0; i<current.nRoots();i++)
//	  {		  
////		  countVec[i]=MathUtils.nint(PathCount.logPathCount(current.getSubtree(i).topology()));
//		  result+=PathCount.logPathCount(current.getSubtree(i).topology());
//		  total+=countVec[i];
////		  result-=MathUtils.logFactorial((int) countVec[i]);
//	  }
////	  result+=MathUtils.logFactorial((int) total);
//	return result;  
//  }
  
  public static long logFact(int n)
  {
    long answer = 1;    // Accumulate the product

    for (long i = 1; i <= n; i++)
      answer *= i;
    
    return answer;
  }

  
  public static double logCountRatio(PartialCoalescentState current, int i0, int i1)
  {
	  int nRoots=current.nRoots();
	  List<Arbre<CoalescentNode>> roots= current.getRoots();
	  double denom=0;
	  for(int i=0;i<nRoots;i++)
	  {
//		  System.out.println(roots.get(i).nLeaves()-1);
		  denom+=roots.get(i).nLeaves()-1;
	  }
	  double numerator=roots.get(i0).nLeaves()+roots.get(i1).nLeaves()-1;
	  if(denom==0)return 0.0;
	  return Math.log(numerator)-Math.log(denom+1);
  }
  public static double nChoose2(double n) { return n*(n-1)/2; }

}

