package pty.learn;
import java.io.*;
import java.util.*;

import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

import pty.Train;
import pty.smc.PartialCoalescentState;
import pty.smc.models.CTMC;
import pty.smc.models.CachedEigenDecomp;

import fig.basic.NumUtils;
import fig.basic.Pair;
import fig.prob.SampleUtils;
import goblin.Taxon;

import nuts.io.IO;
import nuts.math.GMFct;
import nuts.math.MtxUtils;
import nuts.math.RateMtxUtils;
import nuts.math.Sampling;
import nuts.tui.Table;
import nuts.util.Arbre;
import nuts.util.MathUtils;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

import ma.RateMatrixLoader;

/**
 * The following encoding is used in this class for suff stats A
 * 
 * for (a != b):
 * A[a][b] = N(a->b)
 * for (a == b):
 * A[a][a] = T(a)
 * 
 * N(a->b) is the number of transitions from a to b
 * T(a) is the time spent at time in state a 
 * 
 * @author bouchard
 *
 */
public class CTMCExpectations
{
  
  /**
   * for (a != b):
   * A[i][j][a][b] = E[N(a->b)|X_0=i, X_T=j]
   * for (a == b):
   * A[i][j][a][a] = E[T(a)|X_0=i, X_T=j]
   * where 
   * N(a->b) is the number of transitions from a to b in the interval [0,T] and
   * T(a) is the time spent at time in state a in the interval [0,T]
   * @param T
   * @param ed
   * @return suff stat mtx A
   */
  public static double[][][][] expectations(double T, final CachedEigenDecomp ed)
  {
    
    double [] imag = ed.getImagEigenvalues();
    for (int i = 0; i < imag.length; i++)
      if (imag[i] != 0.0)
        throw new RuntimeException();
    
    final int size = ed.getV().getRowDimension();
    // P_T = exp(T Q)
    final double [][] PT = MtxUtils.exp(ed.getV(), ed.getVinv(), ed.getD().times(T)).getArray();
      //MtxUtils.exp(ed.getV(), ed.getD().times(T)).getArray();
    
    double [][] J = computeJ(T, ed);
    double [][] u = ed.getV().getArray();
    double [][] ui = ed.getV().inverse().getArray();
    final double [][] q = ed.getV().times(ed.getD()).times(new Matrix(ui)).getArray();
    
    double [][][][] result = new double[size][size][size][size];
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        for (int a = 0; a < size; a++)
          for (int b = 0; b < size; b++)
          {
            if (a != b)
            {
              final double analyticConvolution = convolution(T, u,ui,   i,a,b,j    ,J);
              
//              {
//                IO.warnOnce("Slow test running");
//                final double numericalConvolution = testConvolution(T, ed, i,a,b,j);
//                if (!MathUtils.close(analyticConvolution, numericalConvolution))
//                  System.out.println("WARNING 1");
//              }
              
              double current = q[a][b] *  analyticConvolution / PT[i][j];
              result[i][j][a][b] = current;
              if (MathUtils.close(0, current) && current < 0)
                current = 0;
              if (current < 0)
                throw new RuntimeException();
            }
            else
            {
              final double analyticConvolution = convolution(T, u,ui,   i,a,a,j    ,J);
              
//              {
//                IO.warnOnce("Slow test running");
//                final double numericalConvolution = testConvolution(T, ed, i,a,a,j);
//                if (!MathUtils.close(analyticConvolution, numericalConvolution))
//                  System.out.println("WARNING 2");
//              }
              
              double current = analyticConvolution   / PT[i][j];
              if (MathUtils.close(T, current) && current > T)
                current = T;
              if (MathUtils.close(0, current) && current < 0)
                current = 0;
              if (current < 0 || current > T)
                throw new RuntimeException();
              result[i][j][a][a] = current;
            }
          }
    return result;
  }
  public static double[][][][] expectations(double T, double [][] rate)
  {
    return expectations(T, new CachedEigenDecomp(new Matrix(rate).eig()));
  }
  
  private static double testConvolution(final double T, final CachedEigenDecomp ed, final int a, final int b, final int c, final int d)
  {
    double result1 = -1;
    {
      UnivariateRealFunction f = new UnivariateRealFunction() {
        
        @Override
        public double value(double t) throws FunctionEvaluationException
        {
          return MtxUtils.exp(ed.getV(), ed.getVinv(), ed.getD().times(t))  .getArray()[a][b] *
                 MtxUtils.exp(ed.getV(), ed.getVinv(), ed.getD().times(T-t)).getArray()[c][d];
        }
      };
      UnivariateRealIntegrator integrator = new TrapezoidIntegrator(f);
      
      try
      {
        result1= integrator.integrate(0,T);
      } catch (Exception e)
      {
        throw new RuntimeException();
      }
    }
    double result2 = -1;
    
    // try again more expensive but one step closer to analytic...
    {
      final double [][] u = ed.getV().getArray();
      final double [][] ui = ed.getV().inverse().getArray();
      final double [] ev = ed.getRealEigenvalues();
      
      UnivariateRealFunction f = new UnivariateRealFunction() {
        
        @Override
        public double value(double t) throws FunctionEvaluationException
        {
          double sum1 = 0.0;
          for (int i = 0; i< u.length; i++)
            sum1 += Math.exp(t * ev[i]) * u[a][i] * ui[i][b];
          
          {
            double otherMethod = MtxUtils.exp(ed.getV(), ed.getVinv(), ed.getD().times(t))  .getArray()[a][b];
            
            if (!MathUtils.close(otherMethod, sum1))
              System.out.println("Divergence:" + otherMethod + "," + sum1);
          }
          
          {
            
          }
          
          double sum2 = 0.0;
          for (int j = 0; j < u.length; j++)
            sum2 += Math.exp((T-t)*ev[j]) * u[c][j] * ui[j][d];
          
          {
            double otherMethod2 = MtxUtils.exp(ed.getV(), ed.getVinv(), ed.getD().times(T-t)).getArray()[c][d];
            
            if (!MathUtils.close(otherMethod2, sum2))
              System.out.println("Divergence2:" + otherMethod2 + "," + sum2);
          }
          
          return sum1 * sum2;
        }
      };
      UnivariateRealIntegrator integrator = new TrapezoidIntegrator(f);
      try
      {
        result2= integrator.integrate(0,T);
      } catch (Exception e)
      {
        throw new RuntimeException();
      } 
    }
    
    if (!MathUtils.close(result1, result2))
      System.out.println("Divergence in result" + result1+"," + result2);

    
    return result1;
  }
  
  /**
   * computes:
   * \int_0^T P_{a->b}(t) P_{c->d}(T-t) dt,
   * where 
   * P(t) = exp(t Q)
   * TODO: might be able to express as matrix multiplication if bottleneck
   */
  private static double convolution(double T, double [][] u, double[][] ui,
      int a, int b, int c, int d, double [][] J)
  {
    double sum = 0.0;
    final int size = u.length;
    for (int i = 0; i < size; i++)
    {
      double subSum = 0.0;
      for (int j = 0; j < size; j++)
        subSum += u[c][j]*ui[j][d]*J[i][j];
      sum += subSum * u[a][i] * ui[i][b];
    }
    return sum;
  }
  

  
  /**
   * computes: e^{T \lambda_j} \int_0^T exp(t(\lambda_i-\lambda_j)) dt
   * @param T
   * @param ed
   * @return
   */
  private static double [][] computeJ(double T, CachedEigenDecomp ed)
  {
    for (double x : ed.getImagEigenvalues())
      if (x != 0.0)
        throw new RuntimeException();
    final double [] ev = ed.getRealEigenvalues();
    final int size = ev.length;
    double [][] result = new double[size][size];
    for (int j = 0; j < size; j++)
    {
      final double factor = Math.exp(T*ev[j]);
      if (factor > 0.0)
        for (int i = 0; i < size; i++)
        {
          final double dLambda = ev[i] - ev[j];
          if (MathUtils.close(dLambda, 0.0))
            result[i][j] = factor * T;
          else
            result[i][j] = factor * (Math.exp(T*dLambda)-1.0) / dLambda;
        }
    }
    return result;
  }
  
  /*
   * Below this is a set of test cases that simulate CTMCs and check
   * the above integrals are correct
   */

  public static Matrix getSufficientStatistics(List<Pair<Integer,Double>> traj, int nChars)
  {
    Matrix result = new Matrix(nChars,nChars);
    int prev = -1;
    for (Pair<Integer,Double> seg : traj)
    {
      int cur = seg.getFirst();
      result.getArray()[cur][cur] += seg.getSecond();
      if (prev != -1)
        result.getArray()[prev][cur]++;
      prev = cur;
    }
    return result;
  }
  public static List<Pair<Integer,Double>> simulate(
      double T, Random rand, double [][] Q)
  {
    double [] sd = RateMtxUtils.getStationaryDistribution(Q);
    int initState = SampleUtils.sampleMultinomial(rand, sd);
    return simulate(initState, T, rand, Q);
  }
  public static List<Pair<Integer,Double>> simulate(int startState,
      double T, Random rand, double [][] Q)
  {
    double [][] embeddedMarkovChain = RateMtxUtils.getJumpProcess(Q);
    List<Pair<Integer,Double>> result = new ArrayList<Pair<Integer,Double>>();
    double totalTime = 0.0;
    int state = startState;
    while (totalTime < T)
    {
      // sample waiting time
      double currentRate = -Q[state][state];
      if (currentRate == 0) // absorbing state !
      {
        result.add(Pair.makePair(state, Double.POSITIVE_INFINITY));
        return result;
      }
      double time = Sampling.sampleExponential(rand, 1.0/currentRate);
      
      totalTime += time;
      if (totalTime > T) time = time - (totalTime -T);
      result.add(Pair.makePair(state, time));
      // transition!
      state = SampleUtils.sampleMultinomial(rand, embeddedMarkovChain[state]);
    }
    return result;
  }
  public static List<List<Pair<Integer,Double>>> simulate(int startState, int endState, int nTrials,
      double minTime, Random rand, double [][] Q)
  {
    List<List<Pair<Integer,Double>>> result = new ArrayList<List<Pair<Integer,Double>>>();
    for ( int i = 0; i < nTrials; i++)
    {
      List<Pair<Integer,Double>> current = simulate(startState, minTime, rand, Q);
      if (current.get(current.size()-1).getFirst() == endState)
        result.add(current);
//      System.out.println(toString(current,1)+"\n---\n");
      
    }
    return result;
  }
  public static int stateAtT(List<Pair<Integer,Double>> seq, double time)
  {
    double current = 0;
    for (Pair<Integer,Double> seg : seq)
    {
      int state = seg.getFirst();
      current += seg.getSecond();
      if (time <= current + 0.00001)
        return state;
    }
    return -1;
  }
  public static double visitTime(List<Pair<Integer,Double>> seq, int state)
  {
    double sum = 0.0;
    for (Pair<Integer,Double> seg : seq)
      if (seg.getFirst() == state)
        sum += seg.getSecond();
    return sum;
  }
  public static int nTransitions(List<Pair<Integer,Double>> seq, int a, int b)
  {
    int sum = 0;
    int prev = seq.get(0).getFirst();
    for (int i = 1; i < seq.size(); i++)
    {
      int current = seq.get(i).getFirst();
      if (a == prev && b == current)
        sum++;
      prev = current;
    }
    return sum;
  }
  public static String toString(List<Pair<Integer,Double>> seq, double timeIncr)
  {
    StringBuilder result = new StringBuilder();
    int cState; double cTime = 0;
    while ((cState = stateAtT(seq, cTime)) != -1)
    {
      for (int i = 0; i < cState; i++)
        result.append(" ");
      cTime += timeIncr;
      result.append("|\n");
    }
    return result.toString();
  }
  
  public static void testEM()
  {
    double bl = 0.1;
    Random rand = new Random(1);
    final int nObservations = 1000,
              nEmIters = 100;
    double [] gen_sd = new double[]{0.4,0.6},
              init_sd= new double[]{0.4,0.6};
    double gen_rate = 0.1, init_rate = 0.1;
    final int nChars = gen_sd.length;
    double [][] gen_q = RateMtxUtils.reversibleRateMtx(gen_rate, gen_sd),
                init_q= RateMtxUtils.reversibleRateMtx(init_rate,init_sd);
    
    System.out.println("Gen:\n" + Table.toString(gen_q));
    System.out.println("sd: " + Arrays.toString(RateMtxUtils.getStationaryDistribution(gen_q) )+ "\n");
    System.out.println("Init:\n" + Table.toString(init_q));
    System.out.println("sd: " + Arrays.toString(RateMtxUtils.getStationaryDistribution(init_q) )+ "\n");
//    CTMC gen = new CTMC.SimpleCTMC(WalsExperiment.binaryRateMtx(1.0,gen_sd), 1);
    
//         init= new CTMC
    // generate observations from the generating distribution: length, pair of states
    List<Pair<Double,Pair<Integer,Integer>>> observations = generateObservations(gen_q, nObservations, rand,bl);
    
    double [][] currentParams = init_q;
    for (int emIter = 0; emIter < nEmIters; emIter++)
    {
      double [][] p = RateMtxUtils.marginalTransitionMtx(currentParams, bl);
      double [] stat = RateMtxUtils.getStationaryDistribution(currentParams);
      // estimate expected suff stats
      CachedEigenDecomp ed = new CachedEigenDecomp(new Matrix(currentParams).eig());
      Matrix suffStats = new Matrix(nChars,nChars);
      for (Pair<Double,Pair<Integer,Integer>> obs : observations)
      {
        double [][][][] expss = expectations(obs.getFirst(),ed);
        int s1 = obs.getSecond().getFirst();
        int s2 = obs.getSecond().getSecond();
        suffStats.plusEquals(new Matrix(expss[s1][s2]));
//        double [] post = new double[2];
//        for (int root = 0; root < nChars; root++)
//          post[root] = stat[root]*p[root][observedLeaf];
//        NumUtils.normalize(post);
//        for (int root = 0; root < nChars; root++)
//          suffStats.plusEquals(new Matrix(expss[root][observedLeaf]).times(post[root]));
      }
      System.out.println("Expected suffstats:\n" + Table.toString(suffStats));
      // estimate new params
      currentParams = Estimators.getGeneralRateMatrixMLE(suffStats);
      System.out.println("After iter " + (emIter+1) + "/" + nEmIters + ":\n" + Table.toString(currentParams));
      System.out.println("sd: " + Arrays.toString(RateMtxUtils.getStationaryDistribution(currentParams) )+ "\n");
    }
    
  }
  
  private static List<Pair<Double, Pair<Integer, Integer>>> generateObservations(
      double[][] Q, int observations, Random rand, double bl)
  {
    List<Pair<Double,Pair<Integer,Integer>>> result = new ArrayList<Pair<Double,Pair<Integer,Integer>>>();
    final int nChars = Q.length;
    Matrix sampledSuffStats = new Matrix(nChars,nChars);
    double [] sd = RateMtxUtils.getStationaryDistribution(Q);
    int prevState = SampleUtils.sampleMultinomial(rand, sd);
//    double visitAt0 = 0.0;
    for (int i =0; i < observations; i++)
    {
      double T = bl; //2*rand.nextDouble();
      List<Pair<Integer,Double>> traj = simulate(prevState, T, rand, Q);
//      visitAt0 += visitTime(traj, 0);
      sampledSuffStats.plusEquals(getSufficientStatistics(traj, nChars));
      int initState = traj.get(0).getFirst(),
          lastState = traj.get(traj.size()-1).getFirst();
      result.add(Pair.makePair(T, Pair.makePair(initState,lastState)));
      prevState = lastState;
    }
    System.out.println("Sampled suffstats:\n" + Table.toString(sampledSuffStats)); 
//    System.out.println("Visit at zero:" + visitAt0);
    return result;
  }
  

  
  public static void main(String [] args)
  {
    
//    for (int s = 2; s < 100; s++)
//    {
//      double [] sd = new double[s];
//      for (int i = 0; i < s; i++)
//        sd[i] = 1.0 / ((double) s);
//      double [][] rateMtx = RateMtxUtils.reversibleRateMtx(1.0,sd);
//      long start = System.currentTimeMillis();
//      expectations(1.0, rateMtx);
//      System.out.println("" + s + "\t" + ((double) (System.currentTimeMillis() - start))/ 1000.0 + "s");
//    }
    
    
    
    testEM();
    if (true) return;
//    if (true) return;
//    double [][] prot = RateMatrixLoader.hky85(); //dayhoff();
//    for (int i = 0; i < 100*100*1000; i++)
//    {
//      double [][][][] test = expectations(1.0, new Matrix(prot).eig());
//      if (i %10000 == 0 )
//      System.out.print(".");
//    }
//    
//    System.out.println("\n---");
    ///
    Random rand = new Random(1);
    double [] sd = new double[] { 0.5, 0.5};
    double rate = 1.0;
    double [][] rateMtx = RateMtxUtils.reversibleRateMtx(rate,sd); // // RateMatrixLoader.hky85(); // RateMatrixLoader.dayhoff();//
    sd = RateMtxUtils.getStationaryDistribution(rateMtx);
    int nTrials = 1000;
    double T = 10;
    
//    RateMtxUtils.checkReversibleRateMtx(rateMtx);
//    //
//    // three ways of computing stat: 
//    // 2 bogus, 1 correct!
//    double [] sd1 = MtxUtils.topEigenvector(RateMatrixLoader.marginalTransitionMtx(dna,1));
//    System.out.println(Arrays.toString(sd1));
//    double [] sd2 = MtxUtils.topEigenvector(dna);
//    System.out.println(Arrays.toString(sd2));
    System.out.println(Arrays.toString(RateMtxUtils.getStationaryDistribution(rateMtx)));
//    //
    
    Matrix ss = new Matrix(rateMtx.length,rateMtx.length);
    for (int i = 0; i < nTrials; i++)
    {
      List<Pair<Integer,Double>> traj = simulate(T, rand, rateMtx);
      ss.plusEquals(getSufficientStatistics(traj, rateMtx.length));
    }
    System.out.println("Gold:\n" + Table.toString(rateMtx));
    System.out.println("Est:\n" + Table.toString(Estimators.getGeneralRateMatrixMLE(ss)));
    
    System.out.println("---");
    int startState = 0; int endState = 0, visit = 0;

    int a = 0, b = 1;
    
    List<List<Pair<Integer,Double>>> samples =
      simulate(startState, endState, nTrials, T, rand, rateMtx);
    double mean = 0.0;
    double meanTr = 0.0;
    
    for (List<Pair<Integer,Double>> sample : samples)
    {
//      System.out.println(toString(sample, 1));
      mean += visitTime(sample, visit);
      meanTr += nTransitions(sample,a,b);
    }
    System.out.println("Simulation visit=" + (mean/((double) samples.size())));
    
    double [][][][] analytic = expectations(T, new CachedEigenDecomp(new Matrix(rateMtx).eig()));
//    System.out.println(Arrays.deepToString(analytic));
    System.out.println("Analytic visit=" + analytic[startState][endState][visit][visit]);
    
    System.out.println("---");
    
    System.out.println("Simulation tr=" + (meanTr/((double) samples.size())));
    
    System.out.println("Analytic tr=" + analytic[startState][endState][a][b]);
    
//    double nTrials = 100000;
//    double nFound = 0.0;
//    Random rand = new Random(1);
//    for (int i = 0; i<nTrials; i++)
//    {
//      double [][] d= new double[4][4];
//      for (int r = 0; r < 4; r++)
//      {
//        double sum = 0.0;
//        for (int c = 0; c < 4; c++)
//          if (c != r)
//          {
//            double n = rand.nextDouble();
//            d[r][c] = n;
//            sum += n;
//          }
//        d[r][r] = -sum;
//      }
//      EigenvalueDecomposition ed = new Matrix(d).eig();
//      searchEv:for (int c = 0; c < 4; c++)
//        if (!MathUtils.close(0.0, ed.getImagEigenvalues()[c]))
//        {
//          System.out.println("Found:\n" + Table.toString(new Matrix(d)) + "has im ev" + Arrays.toString(ed.getImagEigenvalues()));
//          nFound++;
//          break searchEv;
//        }
//    }
//    System.out.println("Density:" + (nFound/nTrials));
  }
      
//  public static void main(String [] args)
//  {
//    double [][] dayhoff = RateMatrixLoader.dayhoff();
//    Matrix m = new Matrix(dayhoff);
//    EigenvalueDecomposition ed = m.eig();
//    System.out.println("im="+Arrays.toString(ed.getImagEigenvalues()));
//    // check  exp Q = exp ( U D U^T) = U exp(D) U^T
////    System.out.println(Table.toString(dayhoff));
//    Matrix test = ed.getV().times(ed.getD()).times(ed.getV().transpose());
////    System.out.println(Table.toString(test));
//    System.out.println(Table.toString(MtxUtils.exp(ed.getV(),ed.getD())));
//    System.out.println(Table.toString(
//        ed.getV().times(new Matrix(MtxUtils.exp(ed.getD().getArray()))).times(ed.getV().transpose())
//        ));
//  }
  
  
}
