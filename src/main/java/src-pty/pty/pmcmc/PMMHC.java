package pty.pmcmc;


import java.util.*;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.io.TreeEvaluator.TreeMetric;
import pty.smc.NCPriorPriorKernel;

import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.PriorPriorKernel;
import pty.smc.models.CTMC;

import gep.util.OutputManager;
import ev.poi.processors.TreeDistancesProcessor;
import nuts.math.Sampling;
import nuts.util.CollUtils;


public class PMMHC 
{
	private final Dataset dataset;
	private final ParticleFilter<PartialCoalescentState> pf;
	private ParticleKernel<PartialCoalescentState> kernel=null;
	private final TreeDistancesProcessor tdp;
	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
	private RootedTree currentSample = null;
	private double trans2tranv;
	public double a=1.5;
	public static OutputManager outMan = new OutputManager();
	//	private double nonclockTreeRate;	
	public boolean stop=false;
	private int tryCountBeforAccept=0;
//	double[] metricResults=new double[TreeEvaluator.coreTreeMetrics.size()];
//	public  UnrootedTree goldut;
	public String out="PMMH";
	private int iter=0;


	public PMMHC(Dataset dataset0, ParticleFilter<PartialCoalescentState> pf,  TreeDistancesProcessor tdp, double trans2tranv0)
	{
		this.trans2tranv=trans2tranv0;
		this.dataset=dataset0;
		this.pf = pf;
		this.tdp = tdp;
		//		this.nonclockTreeRate=nonclockTreeRate;
	}

	public void next(Random rand)
	{
		iter++;
		//		System.out.println(rand.nextDouble());
		RootedTree previousSample = currentSample;		
		double scale=Sampling.nextDouble(rand, 1.0/a, a);
		//		System.out.println(scale);
		double proposedTrans2tranv = scale*trans2tranv; 
		//		double proposedTrans2tranv = 2.0;
		// sample from PF
		StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();		 
//		TreeDistancesProcessor innerTdp = new TreeDistancesProcessor();
////		List<ParticleProcessor<PartialCoalescentState>> processors =CollUtils. list();
//		processors.add(pro);
//		processors.add(innerTdp);
		double acceptPr = 0.0;
//		double[] currentMetricResults=new double[TreeEvaluator.coreTreeMetrics.size()];
//		double[] metricVec=new double[TreeEvaluator.coreTreeMetrics.size()];		
		try 
		{
			tryCountBeforAccept++;
//			if(tryCountBeforAccept>10) stop=true;

			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);  
			PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, ctmc, true);  // clock	          					
			kernel = new PriorPriorKernel(init);
			//			kernel.nonclockTreeRate=nonclockTreeRate;			
//			if(tryCountBeforAccept>5 && pf.N<10000)
//			{
//				pf.N=pf.N+500;
//				tryCountBeforAccept=0;
//			}
//			pf.sample(kernel, processors);
			pf.sample(kernel, pro);
			// compute accept/reject
			double marginalLoglike=pf.estimateNormalizer()+init.logLikelihood();
			final double logRatio =marginalLoglike  - previousLogLLEstimate + 0.5*trans2tranv*(1-scale) + Math.log(scale);
//						System.out.println("log(scale) "+Math.log(scale)+" proposed MarginalLike "+pf.estimateNormalizer()+" Currentlog MarginalLike: "+previousLogLLEstimate+"; log Ratio: "+logRatio);
			acceptPr = Math.min(1, Math.exp(logRatio));
//			System.out.println(pf.N+" "+acceptPr+": "+marginalLoglike  +" -"+ previousLogLLEstimate+(0.5*trans2tranv*(1-scale)+Math.log(scale)));
			if (currentSample != null &&  Double.isInfinite(acceptPr))
				throw new RuntimeException();
//			UnrootedTree  inferred = innerTdp.getConsensus();				
//			for (int i=0; i<TreeEvaluator.coreTreeMetrics.size(); i++)
//			{
//				currentMetricResults[i] = TreeEvaluator.coreTreeMetrics.get(i).score(inferred, goldut);
//				metricVec[i]=acceptPr*currentMetricResults[i]+(1-acceptPr)*metricResults[i];
//			}

			final boolean accept = Sampling.sampleBern(acceptPr, rand);
			if (accept)
			{
//				tryCountBeforAccept=0;
				// sample from sample
				PartialCoalescentState sampled = pro.sample(rand);				
				currentSample = sampled.getFullCoalescentState();
				// set ll!
				previousLogLLEstimate = marginalLoglike;
				trans2tranv=proposedTrans2tranv;
//				for (int i=0; i<TreeEvaluator.coreTreeMetrics.size(); i++) metricResults[i] = currentMetricResults[i];
			}
		}
		catch (Exception e)
		{
			// total fail!
		}

//				System.out.println(currentSample.topology());
		// update tdp
		tdp.process(currentSample);
		// log some stats
		final int tSize = currentSample.topology().nLeaves();
		outMan.write(out,
				"Iter", iter,
				"treeSize", tSize, 
				"pfN",pf.N, 
				"acceptPr", acceptPr, 
				//        "maskSparsity", currentSparsity,
				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
				"trans2tranv", trans2tranv, 
				"LogLikelihood", previousLogLLEstimate
//				TreeEvaluator.coreTreeMetrics.get(0), metricVec[0],
//				TreeEvaluator.coreTreeMetrics.get(1), metricVec[1],
//				TreeEvaluator.coreTreeMetrics.get(2), metricVec[2],
//				TreeEvaluator.coreTreeMetrics.get(3), metricVec[3],
//				TreeEvaluator.coreTreeMetrics.get(4), metricVec[4]
				);
	}
}
