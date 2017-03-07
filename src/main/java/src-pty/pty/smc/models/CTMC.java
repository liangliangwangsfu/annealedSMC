/**
 * 
 */
package pty.smc.models;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import conifer.evol.GTR;
import dr.math.distributions.GammaDistribution;
//import dr.math.distributions.GammaDistribution;
//import math.distributions.GammaDistribution;

import pty.ObservationDimensions;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

import fig.basic.NumUtils;
import fig.prob.Gamma;

import ma.RateMatrixLoader;
import ma.SequenceType;
import nuts.math.MtxUtils;
import nuts.math.RateMtxUtils;
import nuts.tui.Table;
import nuts.util.MathUtils;

/**
 * Continuous Time Markov Chain with a finite set of observations
 * 
 * Allows different state space/params for different sites
 * 
 * @author bouchard
 *
 */
public interface CTMC extends Serializable, ObservationDimensions
{
  
  /**
   * @param site
   * @param t
   * @return top state -> bottom state prs
   */
  public double [][] getTransitionPr(int site, double t);
  public double [] getInitialDistribution(int site);
  public CachedEigenDecomp getRateMtx(int site);
  public boolean isSiteTied();
//  public boolean isReversible();
  
//  public static class RateMatrix
//  {
//    public double [] getStationaryDistribution() { return null; }
//    public double [][] getTransitionPr(double t) { return null; }
//  }
  
  public static final class GeneralCTMC implements CTMC
  {
    private static final long serialVersionUID = 1L;
    private final List<CachedEigenDecomp> Qs = new ArrayList<CachedEigenDecomp>();
    private final List<Matrix> originalQs = new ArrayList<Matrix>();
    private final List<double[]> statDistn = new ArrayList<double[]>();
    public GeneralCTMC(List<double[][]> Qs)
    {
      for (double [][] Q : Qs)
      {
        originalQs.add(new Matrix(Q));
        double [] stat = RateMtxUtils.getStationaryDistribution(Q);
        MathUtils.checkIsProb(stat);
        statDistn.add(stat);
        this.Qs.add(new CachedEigenDecomp(new Matrix(Q).eig()));
      }
    }
    
    
    public double[] getInitialDistribution(int site)
    {
      return statDistn.get(site);
    }
    public double[][] getTransitionPr(int site, double t)
    {
      return RateMtxUtils.marginalTransitionMtx(Qs.get(site).getV(), Qs.get(site).getVinv(), Qs.get(site).getD(), t);
    }
    public int nCharacter(int site)
    {
      return statDistn.get(site).length;
    }
    public int nSites()
    {
      return Qs.size();
    }
    public CachedEigenDecomp getRateMtx(int site)
    {
      return Qs.get(site);
    }
    @Override
    public String toString()
    {
      StringBuilder result = new StringBuilder();
      result.append("GeneralCTMC:\n");
      for (int i =0 ; i < nSites(); i++)
        result.append("Site " + i + ":\n" + Table.toString(originalQs.get(i)) + "\n");
      return result.toString();
    }
    public boolean isSiteTied()
    {
      return false;
    }
  }
  
  /**
   * Each site is assumed to have the same state space and params,
   * initial distn is the stat. distn of the rate mtx
   * @author bouchard
   */
  public static final class SimpleCTMC implements CTMC
  {
    private static final long serialVersionUID = 1L;
    private final CachedEigenDecomp Q;
    private final Matrix originalQ;
    private final double[] statDistn;
    private final int nSites;
    public static SimpleCTMC dnaCTMC(int nSites)
    {
      return new SimpleCTMC(RateMatrixLoader.k2p(), nSites);
    }
    
    public static SimpleCTMC dnaCTMC(int nSites, double trans2tranv)                // added by Liangliang on July 8, 2011
    {
      return new SimpleCTMC(RateMatrixLoader.k2p(trans2tranv), nSites);
    }
    
    
    public static SimpleCTMC proteinCTMC(int nSites)
    {
      return new SimpleCTMC(RateMatrixLoader.dayhoff(), nSites);
    }
    public static SimpleCTMC fromSequenceType(int nSites, SequenceType st,double scale)
    {
//      if      (st == SequenceType.DNA) return dnaCTMC(nSites); this needs some handling of scale+remove the else comment below
//      else if (st == SequenceType.PROTEIN) return proteinCTMC(nSites);
     /*else*/ if (st == SequenceType.BINARY)
      {
        double [][] rate = new double[][] {{-scale,+scale},{+scale,-scale}};
        return new SimpleCTMC(rate, nSites);
      }
      else throw new RuntimeException();
    }
    public SimpleCTMC(double [][] rate, int nSites) 
    {
      this.nSites = nSites;
      this.originalQ = new Matrix(rate);
      this.Q = new CachedEigenDecomp(originalQ.eig());
      this.statDistn = RateMtxUtils.getStationaryDistribution(rate);
      MathUtils.checkIsProb(this.statDistn);
    }
    public double[] getInitialDistribution(int site)
    {
      return statDistn;
    }
    public double[][] getTransitionPr(int site, double t)
    {
      return RateMtxUtils.marginalTransitionMtx(Q.getV(), Q.getVinv(), Q.getD(), t);
    }
    public int nCharacter(int site) { return statDistn.length; }
    public int nSites() { return nSites; }
    public CachedEigenDecomp getRateMtx(int site)
    {
      return Q;
    }
    @Override
    public String toString()
    {
      return "SimpleCTMC:\n" + Table.toString(originalQ);
    }
    public boolean isSiteTied()
    {
      return true;
    }
  }
  
  /**
   * Each site is assumed to have the same state space and params,
   * initial distn is the stat. distn of the rate mtx
   * @author Liangliang Wang
   */
  public static final class GTRIGammaCTMC implements CTMC
  {
    private static final long serialVersionUID = 1L;
    private final CachedEigenDecomp Q;
    private final Matrix originalQ;
    private final double[] statDistn;
    private final int nSites;
    private final int nGammaCat; 
    private final double pInv; 
//     
    /*
     * stat: stationary state frequencies.
     * subRates: rates of substitutions
     * n: 4 for DNA; i.e. the size of {A,C,G,T}
     * alpha: shape parameters in the Gamma distribution for r, where r is in exp(rQb). 
     * nCategories: number of categories in the discrete Gamma. Typically it is chosen to be 4. 
     * pInv: the proportion of invariant sites. 
     * */
    public GTRIGammaCTMC(double [] stat, double [] subRates, int n, int nSites, double alpha, int nGammaCat, double pInv) 
    {
    	  double [] r=calculateCategoryRates(nGammaCat,alpha, pInv);
    	  int dim=nGammaCat*n;
    	  if(pInv>0 && pInv<1) dim=dim+n;
      double [][] rate0=new double[dim][dim];      	
//    	  double [][] rateMtx=GTR.gtrFromOverParam(stat, subRates, n);
      double [][] rateMtx=GTR.scaleGTRrateMat(stat, GTR.gtrFromOverParam(stat, subRates, n));
     
    	  for (int c = 0; c < nGammaCat; c++)    	      		    
    		  for(int col = 0;  col< n; col++)
    	        for (int row = 0; row < n; row++)    	        	
    	            rate0[c*n+row][c*n+col] = r[c]*rateMtx[row][col]; 

    	  
    	  //  	  System.out.println("Q : "+Table.toString(rate));		  
    	  this.statDistn = new double[dim];//RateMtxUtils.getStationaryDistribution(rate);
//    	  this.statDistn = RateMtxUtils.getStationaryDistribution(rate0);
    	  for (int c = 0; c < nGammaCat; c++)    	      		    
    		  for(int i = 0;  i< n; i++)
    			  this.statDistn[c*n+i] = (1-pInv)*stat[i]/nGammaCat;
//    	  System.out.println(Arrays.toString(this.statDistn));
    	  if(pInv>0 && pInv<1) 
    	  {
    		  for(int i = 0;  i< n; i++)
    			  this.statDistn[n*nGammaCat+i] = pInv*stat[i];		  
    	  }
//    rate0=GTR.scaleGTRrateMat(this.statDistn, rate0);  // I may need to change this when pInv in (0,1)      
    	  this.nSites = nSites;
      this.originalQ = new Matrix(rate0);
      this.Q = new CachedEigenDecomp(originalQ.eig());      
      	  
//      MathUtils.checkIsProb(this.statDistn);
      this.nGammaCat=nGammaCat;
      this.pInv=pInv; 
   }
    
    
    public GTRIGammaCTMC(double [] stat, double [] subRates, int n, int nSites,int nGammaCat, double pInv) 
    {
            	
//    	  double [][] rateMtx=GTR.gtrFromOverParam(stat, subRates, n);
      double [][] rateMtx=GTR.scaleGTRrateMat(stat, GTR.gtrFromOverParam(stat, subRates, n));    
//    	double [][] rateMtx=GTR.gtrFromOverParam(stat, subRates, n);
    	  
  //  	  System.out.println("Q : "+Table.toString(rate));		  
    	  //this.statDistn = new double[dim];//RateMtxUtils.getStationaryDistribution(rate);
    	  this.statDistn = RateMtxUtils.getStationaryDistribution(rateMtx);
//      this.statDistn = stat;
//    	  for (int c = 0; c < nGammaCat; c++)    	      		    
//    		  for(int i = 0;  i< n; i++)
//    			  this.statDistn[c*n+i] = (1-pInv)*stat[i]/nGammaCat;
//    	  if(pInv>0 && pInv<1) 
//    	  {
//    		  for(int i = 0;  i< n; i++)
//    			  this.statDistn[n*nGammaCat+i] = pInv*stat[i];		  
//    	  }
//    	  double [][] rate=GTR.scaleGTRrateMat(this.statDistn, rate0);    	  
    	  this.nSites = nSites;
      this.originalQ = new Matrix(rateMtx);
      this.Q = new CachedEigenDecomp(originalQ.eig());            	  
      MathUtils.checkIsProb(this.statDistn);
      this.nGammaCat=nGammaCat;
      this.pInv=pInv; 
   }
    
    /**
     * discretization of gamma distribution with equal proportions in each
     * category
     * @throws IllegalArgumentException 
     * @throws FunctionEvaluationException 
     * @throws MaxIterationsExceededException 
     */
    public static double[] calculateCategoryRates(int nGammaCat, double alpha, double pInv)  {
    	    	                 	
        double[] categoryRates=new double[nGammaCat];
        if(nGammaCat==1)
        {
        	categoryRates[0]=1;
        	return categoryRates;  
        }

//        if(nGammaCat==5)
//        {
//        	categoryRates=new double[]{0.02121, 0.15549, 0.49708, 1.10712, 3.24910};
//        	return categoryRates;  
//        }
        
        double  propVariable = 1.0 - pInv;

//            final double alpha = shapeParameter.getParameterValue(0);
//          double mean = 0.0;
//           
//            for (int i = 0; i < nGammaCat; i++) {
//                categoryRates[i] = GammaDistribution.quantile((2.0*i + 1.0)/(2.0*nGammaCat), alpha, 1.0/alpha);
//            mean += categoryRates[i];
//            }
//
//           mean = (propVariable*mean)/nGammaCat;
//
//            for (int i = 0; i < nGammaCat; i++) {
//                categoryRates[i] /= mean;
//            }
                        
            double[] quantiles=new double[nGammaCat+1];
        	    quantiles[0]=0;
            for (int i = 1; i <= nGammaCat; i++) 
            	quantiles[i] = GammaDistribution.quantile((2.0*i)/(2.0*nGammaCat), alpha, 1.0/alpha);
        	    //quantiles[nGammaCat]=Double.POSITIVE_INFINITY;
            quantiles[nGammaCat]=200;
            
        	    Rmean rmean= new Rmean(alpha); 
            TrapezoidIntegrator trInt=new TrapezoidIntegrator();
            double sum=0;
            for (int i = 0; i < nGammaCat; i++) {           	                          
            try {
//            	System.out.println("quantile "+i+": "+quantiles[i]);
				categoryRates[i] =trInt.integrate(rmean, quantiles[i], quantiles[i+1])*nGammaCat;
				sum+=categoryRates[i];
			} catch (MaxIterationsExceededException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (FunctionEvaluationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
//            System.out.println((i+1)+": "+categoryRates[i]);
            
            }
//            System.out.println("sum is "+sum);
            
            
//            double sumR=0;
//            for (int i = 0; i < nGammaCat; i++) {
//            System.out.println(categoryRates[i]);
//            sumR=sumR+categoryRates[i];
//            }
//            System.out.println(sumR);
            

//        if (muParameter != null) { // Moved multiplication by mu to here; it also
//                                   // needed by double[] getCategoryRates() -- previously ignored
//            double mu = muParameter.getParameterValue(0);
//             for (int i=0; i < categoryCount; i++)
//                categoryRates[i] *= mu;
//        }

//        ratesKnown = true;
        return categoryRates; 
    }
    
    
//   public int getnCategories()
//   {
//	   return nCategories; 
//   }
//    
//   
//   public double getpInv()
//   {
//	   return pInv; 
//   }
//   
   
    public double[] getInitialDistribution(int site)
    {
      return statDistn;
    }
    public double[][] getTransitionPr(int site, double t)
    {    	     	
    	int m=this.statDistn.length;
    	double[][] result=new double[m][m];
		DoubleMatrix A=DoubleMatrix.zeros(m,m);
		for(int i=0;i<m;i++)
			for(int j=0;j<m;j++)		
		    A.put(i,j,  this.originalQ.get(i, j)*t);

		DoubleMatrix expA=MatrixFunctions.expm(A);
		
		for(int i=0;i<m;i++)
			for(int j=0;j<m;j++)		
              result[i][j]=expA.get(i, j); 		
		return result;
//      return RateMtxUtils.marginalTransitionMtx(Q.getV(), Q.getVinv(), Q.getD(), t);
    }
    public int nCharacter(int site) { return statDistn.length; }
    public int nSites() { return nSites; }
    public CachedEigenDecomp getRateMtx(int site)
    {
      return Q;
    }
    

    @Override
    public String toString()
    {
      return "GTRIGammaCTMC:\n" + Table.toString(originalQ);
    }
    public boolean isSiteTied()
    {
      return true;
    }
    
    public static void main(String [] args)
    {
      double []stat = new double[]{0.3, 0.2, 0.2, 0.3};
      
      double [] rates = new double[] {0.26, 0.18, 0.17, 0.15, 0.11, 0.13};
      
//      double [][] rateMtx = gtrFromOverParam(stat, rates, 4);
      int nCategories=4;
      int n=4;     
//	  int dim=nCategories*n;
//      double [][] rate=new double[dim][dim];      	
//	  double [][] rateMtx=GTR.gtrFromOverParam(stat, rates, n);
//	  double [] r=new double[]{10,8,6,7}; 
//	  for (int c = 0; c < nCategories; c++)    	      		    
//		  for(int col = 0;  col< n; col++)
//	        for (int row = 0; row < n; row++)    	        	
//	            rate[c*n+row][c*n+col] = r[c]*rateMtx[row][col];
//
//      
//      double [] stat2 = RateMtxUtils.getStationaryDistribution(rate);
//      
//      System.out.println(Arrays.toString(stat2));
//      
//      System.out.println(Table.toString(rate));
//      
//      
//      System.out.println(Table.toString(RateMtxUtils.marginalTransitionMtx(rate, 100)));
//      
//      
//      System.out.println(MatrixFunctions.expm(new DoubleMatrix(rate).mul(100.0)));
//      
      CTMC ctmc = new GTRIGammaCTMC(stat, rates, 4, 100, 0.25, nCategories, 0);
      double[][] tranMat=ctmc.getTransitionPr(0, 1);
      System.out.println(Table.toString(tranMat));
      //double[]  statDistn=RateMtxUtils.getStationaryDistribution(tranMat);
      
      //double [][] marginal = RateMtxUtils.marginalTransitionMtx(ctmc.getRateMtx(0), 1.0);
      
      
      System.out.println(Arrays.toString(ctmc.getInitialDistribution(0)));
      
      
      
//      double [][] simple = RateMatrixLoader.k2p();
//      double [] stat = RateMtxUtils.getStationaryDistribution(simple);
//      System.out.println(Arrays.toString(stat));
    }
  }

  public class Rmean implements UnivariateRealFunction
  {
    private double alpha=1.0; 
    public Rmean(double alpha){
    	this.alpha=alpha;
    }    
	@Override
	public double value(double r) throws FunctionEvaluationException {		
//		return r*GammaDistribution.pdf(r, alpha, 1.0/alpha);
//        return r*Math.exp(Gamma.logProb(r, alpha, alpha)); 
//		System.out.println("version 1 "+Math.exp(Gamma.logProb(r, alpha, alpha)));
//		System.out.println("version 2 "+GammaDistribution.pdf(r, alpha, 1.0/alpha));
		return r*GammaDistribution.pdf(r, alpha, 1.0/alpha);
	}
	  
  }
  
}