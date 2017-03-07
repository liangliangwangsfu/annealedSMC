/**
 * 
 */
package pty.smc.models;
import java.io.File;
import static nuts.util.CollUtils.list;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.lucene.util.OpenBitSet;

import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Pair;
import fig.prob.SampleUtils;
import goblin.Taxon;
import ma.MSAPoset;
import ma.RateMatrixLoader;
import ma.SequenceType;
import ma.MSAPoset.Column;
import nuts.io.IO;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils;
import nuts.util.CounterMap;
import nuts.util.Indexer;
import nuts.util.MathUtils;
import pepper.Encodings;
import pty.io.Dataset;
import pty.smc.PartialCoalescentState;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.models.CTMC.SimpleCTMC;

/**
 * Likelihood models based on finite state CTMC
 * Has methods to initialize a likelihood model based on observations,
 * to update the likelihoods for a coalescent event.
 * 
 * Uses scalings instead of logs
 * @author bouchard
 */
public final class GTRIGammaFastDiscreteModelCalculator implements LikelihoodModelCalculator, Serializable 
{
	private static final long serialVersionUID = 1L;

	public LikelihoodModelCalculator combine(
			LikelihoodModelCalculator node1, LikelihoodModelCalculator node2, 
			double v1, double v2, boolean avoidBuildCache)
	{
		return (LikelihoodModelCalculator) calculate(node1, node2, v1, v2, false, avoidBuildCache);
	}

	public static GTRIGammaFastDiscreteModelCalculator observation(CTMC ctmc, double [][] initCache, double alpha, int nGammaCat)
	{    
		
		return new GTRIGammaFastDiscreteModelCalculator(initCache, ctmc, alpha, nGammaCat,false);
	}


	public static GTRIGammaFastDiscreteModelCalculator observation(CTMC ctmc, double [][] initCache, double alpha, int nGammaCat, boolean resampleRoot)
	{    
		return new GTRIGammaFastDiscreteModelCalculator(initCache, ctmc, alpha, nGammaCat,resampleRoot);
	}

	/*
	 * site -> character (std scale)
	 */
	private final ArrayList<double [][]> cache; 
	private final double logLikelihood;
    public final ArrayList<double []> itemizedLogLikelihoods;
	public final CTMC ctmc;
	//  private final boolean resampleRoot;
	private final OpenBitSet [] bitVectorVersion;

	private final double alpha; 
	private final int nGammaCat; 
	private final double[] rates;


	//  private final DiscreteModelCalculator tester;

	public ArrayList<double[][]> getCache(double rate)
	{
		return cache;
	}

	// called when combining
	private GTRIGammaFastDiscreteModelCalculator(ArrayList<double[][]> cache, CTMC ctmc,
			double logLikelihood, double alpha, int nGammaCat,boolean resampleRoot, ArrayList<double []> itemizedLogLikelihoods)//, DiscreteModelCalculator tester)
	{
		this.bitVectorVersion = null;
		//    this.resampledVersion = null; //resampleToBitVector(cache);
		//    this.tester = tester;
		this.cache = cache;
		this.ctmc = ctmc;
		this.logLikelihood = logLikelihood;
		this.alpha=alpha;
		this.nGammaCat=nGammaCat;
		this.rates=CTMC.GTRIGammaCTMC.calculateCategoryRates(nGammaCat, alpha, 0);
		this.itemizedLogLikelihoods=itemizedLogLikelihoods; 
		//    this.itemizedLogLikelihoods = itemizedLogLikelihoods;
		//    this.resampleRoot = resampleRoot;
	}

	// called when initialized
	private GTRIGammaFastDiscreteModelCalculator(double[][] cache, CTMC ctmc, double alpha, int nGammaCat, boolean resampleRoot) // call this if you are initializing !
	{
		if (resampleRoot)
			throw new RuntimeException("Discontinued feature");    
		//    IO.warnOnce("Using instrumented FDM... change that!");
		//    this.tester = DiscreteModelCalculator.observation(ctmc,cache);
	    ArrayList<double []> itemized=new ArrayList<double []>(); 
	    for(int i=0;i<ctmc.nSites();i++) itemized.add(new double[nGammaCat]);
		this.ctmc = ctmc;
		this.itemizedLogLikelihoods = itemized;
		this.logLikelihood = initLogLikelihood(cache, ctmc, this.itemizedLogLikelihoods,nGammaCat);
		this.cache = initCache(cache, ctmc,nGammaCat);
		this.bitVectorVersion = observationToBitVector(cache); //resampleToBitVector(this.cache);
		this.alpha=alpha;
		this.rates=CTMC.GTRIGammaCTMC.calculateCategoryRates(nGammaCat, alpha, 0);
		this.nGammaCat=nGammaCat;		
		//    this.resampleRoot = resampleRoot;

	}



	private static ArrayList<double[][]> initCache(double[][] cache, CTMC ctmc, int nRepeats)
	{
		//    return cache;
		if (!ctmc.isSiteTied())
			throw new RuntimeException();
		double [] initD = ctmc.getInitialDistribution(0);
		double[][] result=new double[cache.length][];
		for (int i = 0; i < result.length; i++)
		{
			// this handles unobserved or partially observed r.v.s
			result[i] = new double[cache[i].length];
			for (int j = 0; j < result[i].length; j++)
				result[i][j] = initD[j] * cache[i][j];
			// Note: the line below is important, caused a bug without it
			NumUtils.normalize(result[i]);
		}
			
		ArrayList<double [][]> allResult =list();     
		for(int i=0;i<nRepeats;i++)
		{    	 
			allResult.add(result);
		}
		return allResult;
	}


	private static double initLogLikelihood(double[][] cache, CTMC ctmc, ArrayList<double []> itemized, int nGammaCat)
	{
		if (!ctmc.isSiteTied())
			throw new RuntimeException();

		//    final boolean needItemized = itemized != null;


		double [] initD = ctmc.getInitialDistribution(0);
		double result = 0.0;
		for (int s = 0; s < cache.length; s++)
		{
			for (int c = 0; c < initD.length; c++)
			{
				final double cur = cache[s][c];
				if (cur != 0.0 && cur != 1.0)      
					throw new RuntimeException();
				final double term = cur * Math.log(initD[c]);
				//        if (needItemized)
				         itemized.get(s)[0] += term;
				result += term;
			}		
//			itemized.get(s)[0]=itemized.get(s)[0]-Math.log(nGammaCat);
			for(int k=1;k<nGammaCat;k++) itemized.get(s)[k]=itemized.get(s)[0];
//			System.out.println("s "+s+" itemized.get(s)[0] "+itemized.get(s)[0]);
		}
		return result;
	}

	public double extendLogLikelihood(double delta)
	{
		throw new RuntimeException();
	}


	public static double quickPeek(      
			final LikelihoodModelCalculator node1, final LikelihoodModelCalculator node2, 
			final double v1, final double v2)
	{
		System.out.println("Doing quick peek");
		final GTRIGammaFastDiscreteModelCalculator 
		cache1 = (GTRIGammaFastDiscreteModelCalculator) (node1),
		cache2 = (GTRIGammaFastDiscreteModelCalculator) (node2);
		final double [][] pr = cache1.ctmc.getTransitionPr(0,v1+v2);
		final double [] initD = cache2.ctmc.getInitialDistribution(0);
		final double [] logInit = new double[initD.length];
		for (int i = 0; i < logInit.length; i++)
			logInit [i]= Math.log(initD[i]);

		final int ncs = cache1.ctmc.nCharacter(0);

		double sum = 0.0;
		for (int x = 0; x < ncs; x++)
		{
			final double cur = logInit[x];
			final double [] curAr = pr[x];
			for (int xPrime = 0; xPrime < ncs; xPrime++)
			{
				final long n = OpenBitSet.intersectionCount(cache1.bitVectorVersion[x], cache2.bitVectorVersion[xPrime]);
				sum += n * (cur + Math.log(curAr[xPrime]));
			}
		}
		return sum;
	}

	/**
	 * When combining two trees in a coalescent event,
	 * return a fresh LikelihoodModelCalculator that
	 * sits at the given height
	 * @param node1
	 * @param node2
	 * @param height
	 * @return
	 */
	private final Object calculate(
			final LikelihoodModelCalculator node1, final LikelihoodModelCalculator node2, 
			final double v1, final double v2,
			final boolean isPeek, final boolean avoidBuildCache)
	{

		if (isPeek && !avoidBuildCache)
			throw new RuntimeException();

		// resulting backward scores are indexed: site -> character, in std scale
		ArrayList<double [][]> result = list(); 
		for(int c=0;c<nGammaCat;c++)
			result.add((avoidBuildCache ? null :Dataset.DatasetUtils.createObsArray(ctmc)));
		
		final GTRIGammaFastDiscreteModelCalculator 
		cache1 = (GTRIGammaFastDiscreteModelCalculator) (node1),
		cache2 = (GTRIGammaFastDiscreteModelCalculator) (node2);

//		 final boolean needItemized = !avoidBuildCache && (cache1.itemizedLogLikelihoods != null);
//		    if (needItemized != (cache2.itemizedLogLikelihoods != null))
//		      throw new RuntimeException();

		if (isPeek && cache1.bitVectorVersion != null && cache2.bitVectorVersion != null)
			return quickPeek(node1, node2, v1, v2);

		//    DiscreteModelCalculator dmc = (DiscreteModelCalculator) tester.combine(cache1.tester, cache2.tester, v1, v2, avoidBuildCache);
		//    System.out.println("DMC: " + dmc.logLikelihood());

		    final ArrayList<double []> newItemized = list();		    
	    for(int s=0;s<ctmc.nSites();s++)    newItemized.add(new double[nGammaCat]);
		
		final ArrayList<double [][]> pr1array=list(), pr2array=list();  
		for(int c=0;c<nGammaCat;c++) 
		{
			pr1array.add(ctmc.getTransitionPr(0,v1*rates[c])); 
			pr2array.add(ctmc.getTransitionPr(0,v2*rates[c]));
		}

		final double [] initD = ctmc.getInitialDistribution(0);
		final double [] invD = new double[initD.length];

		for (int i = 0; i < invD.length; i++)
		{
			initD[i]=initD[i]; 
			invD[i] = 1.0/initD[i];
		}

		//    {
			//      System.out.println(Arrays.deepToString(pr1));
		//      System.out.println(Arrays.deepToString(pr2));
		//      System.out.println(Arrays.toString(initD));
		//      System.out.println(Arrays.toString(invD));
		//    }


		final int ns = ctmc.nSites();
		final int ncs = ctmc.nCharacter(0);
		double logNorm = 0.0, tempNorm = 1.0;
		for (int s = 0; s < ns; s++)
		{
			for(int c=0;c<nGammaCat;c++)
			{
				// compute the norm of result[s][x]
				double currentNorm = 0.0;
				// optionally, fill an unormalized version of cache[s][x]...

				final double [] c1 = cache1.cache.get(c)[s], c2 = cache2.cache.get(c)[s];
				for (int x = 0; x < ncs; x++)
				{
					final double [] cPr1 = pr1array.get(c)[x], cPr2 = pr2array.get(c)[x];
					double s1 = 0.0, s2 = 0.0;
					for (int xPrime = 0; xPrime < ncs; xPrime++)
					{
						final double cInvD = invD[xPrime];
						final double cc1 = c1[xPrime]; if (cc1 != 0.0) s1 += cPr1[xPrime] * cInvD * cc1;
						final double cc2 = c2[xPrime]; if (cc2 != 0.0) s2 += cPr2[xPrime] * cInvD * cc2;
					}
					final double currentProd = initD[x] * s1 * s2;
					if (!avoidBuildCache) result.get(c)[s][x]= currentProd;
//					System.out.println(""+ c+" "+rates[c]+" "+currentProd);					
					currentNorm += currentProd;
				}
				if (!avoidBuildCache) 
					MathUtils.normalizeAndGetNorm(result.get(c)[s]);
				//      if (needItemized)
				       newItemized.get(s)[c] = Math.log(currentNorm) + cache1.itemizedLogLikelihoods.get(s)[c] + cache2.itemizedLogLikelihoods.get(s)[c];
//				       System.out.println(Math.log(currentNorm) +" "+ cache1.itemizedLogLikelihoods.get(s)[c] +" "+ cache2.itemizedLogLikelihoods.get(s)[c]);
//						tempNorm *= currentNorm;
//						System.out.println(currentNorm);
			}
			logNorm+=SloppyMath.logAdd(newItemized.get(s))-Math.log(nGammaCat); 
//			System.out.println(SloppyMath.logAdd(newItemized.get(s))-Math.log(nGammaCat));

			
//			double abs = Math.abs(tempNorm);
//			if (abs < 1e-100 || abs > 1e+100 || s == ns -1)
//			{
//				//        {
//				//          System.out.println("\t"+ tempNorm);
//				//        }
//
//				logNorm += Math.log(tempNorm);
//				tempNorm = 1.0;
//			}
		}
		// NB: exp(lognorm) is equal to P(y_{left+right})/P(y_left)/P(y_right)
//		final double newLogLikelihood = logNorm + cache1.logLikelihood + cache2.logLikelihood;
//		System.out.println("logNorm "+logNorm+" cache1.logLikelihood "+cache1.logLikelihood+ "cache2.logLikelihood"+cache2.logLikelihood+":"+newLogLikelihood);
		final double newLogLikelihood = logNorm;
		
//         System.out.println("Ours: " + newLogLikelihood+ " "+ logNorm+" "+cache1.logLikelihood+" "+cache2.logLikelihood);
		
		// hence, to get P(y_{left+right}), can take the product (sum the logs)
		// i.e.: incremental telescoping product

		//    System.out.println("Ours: " + newLogLikelihood);
		//    System.out.println();
		//    MathUtils.checkClose(newLogLikelihood, dmc.logLikelihood());

		if (isPeek) return newLogLikelihood;
		else        return new GTRIGammaFastDiscreteModelCalculator(result, ctmc, newLogLikelihood, alpha, nGammaCat, false,newItemized); 
	}
	
	
//	(ArrayList<double[][]> cache, CTMC ctmc,
//			double logLikelihood, double alpha, int nGammaCat,boolean resampleRoot) 

	private static OpenBitSet [] observationToBitVector(double [][] message)
	{
		final int nChars = message[0].length;
		final int nSites = message.length;
		OpenBitSet [] result = new OpenBitSet[nChars];
		for (int i = 0; i < result.length; i++)
			result[i] = new OpenBitSet(nSites);
		for (int site = 0; site < message.length; site++)
		{
			int idx = findOne(message[site]);
			if (idx == -1)
				return null;
			result[idx].fastSet(site);
		}
		return result;
	}


	private static int findOne(double[] ds)
	{
		int result = -1;
		for (int i = 0; i < ds.length; i++)
		{
			final double d = ds[i];
			if (d == 1)
			{
				if (result != -1)
					return -1; // more than one == 1
				result = i;
			}
		}
		return result;
	}

	public double logLikelihood()
	{
		return logLikelihood;
	}

	public boolean isReversible()
	{
		throw new RuntimeException();
	}

	public double peekCoalescedLogLikelihood(LikelihoodModelCalculator node1,
			LikelihoodModelCalculator node2, double delta1, double delta2)
	{
		return (Double) calculate(node1, node2, delta1, delta2, true, true);
	}
}