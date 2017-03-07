package pty.pmcmc;

import java.util.*;

import pty.RootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.smc.LazyNonclockPriorPrior;
import pty.smc.LazyPCS;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.models.CTMC;

import gep.util.OutputManager;
import ev.poi.processors.TreeDistancesProcessor;
import fig.basic.NumUtils;
import fig.prob.SampleUtils;
import nuts.math.Sampling;


public class LAZYPMMHNC 
{
	private final Dataset dataset;
	private final ParticleFilter<LazyPCS> pf;
	private LazyNonclockPriorPrior kernel=null;
//	private NCPriorPriorKernel kernel=null;
//	private ParticleKernel<PartialCoalescentState> kernel=null;
	private final TreeDistancesProcessor tdp;
	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
	private RootedTree currentSample = null;
	private double trans2tranv;
	public double a=1.2;
	public static OutputManager outMan = new OutputManager();
	private double nonclockTreeRate;
	public String treeName;

	
	public LAZYPMMHNC(Dataset dataset0, ParticleFilter<LazyPCS> pf,  TreeDistancesProcessor tdp, double trans2tranv0, double nonclockTreeRate)
	{
		this.trans2tranv=trans2tranv0;
		this.dataset=dataset0;
		this.pf = pf;
		this.tdp = tdp;
		this.nonclockTreeRate=nonclockTreeRate;
	}

	public void next(Random rand)
	{
//		System.out.println(rand.nextDouble());
		RootedTree previousSample = currentSample;
		double scale=Sampling.nextDouble(rand, 1.0/a, a);
//		System.out.println(scale);
		double proposedTrans2tranv = scale*trans2tranv; 
		//		double proposedTrans2tranv = 2.0;
		// sample from PF
		StoreProcessor<LazyPCS> pro = new StoreProcessor<LazyPCS>();
		double acceptPr = 0.0;
		try 
		{
			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);  
			PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, ctmc, false);	          
			kernel = new LazyNonclockPriorPrior(init);		
//			kernel = new NCPriorPriorKernel(init);
//			kernel.nonclockTreeRate=nonclockTreeRate;
			
			pf.sample(kernel, pro);
			// compute accept/reject

			final double logRatio = pf.estimateNormalizer() - previousLogLLEstimate + Math.log(scale);
			//			System.out.println("log(scale) "+Math.log(scale)+" proposed MarginalLike "+pf.estimateNormalizer()+" Currentlog MarginalLike: "+previousLogLLEstimate+"; log Ratio: "+logRatio);

			acceptPr = Math.min(1, Math.exp(logRatio));
			if (currentSample != null &&  Double.isInfinite(acceptPr))
				throw new RuntimeException();
			final boolean accept = Sampling.sampleBern(acceptPr, rand);
			
			if (accept)
			{
				// sample from sample
				LazyPCS sampled = pro.sample(rand);
				currentSample = sampled.getState().getFullCoalescentState();
				//        List<LazyPCS> particleList = pf.getSamples(); 
				//		double[] finalWeights = pf.getLogWeights(); // log weights. 
				//		
				//		NumUtils.expNormalize(finalWeights);
				//		final int sampledIndex = SampleUtils.sampleMultinomial(rand, finalWeights);
				//		  currentSample =particleList.get(sampledIndex).getState().getFullCoalescentState();

				// set ll!
				previousLogLLEstimate = pf.estimateNormalizer();
				trans2tranv=proposedTrans2tranv;
			}
		}
		catch (Exception e)
		{
			// total fail!
		}

//		System.out.println(currentSample.topology());
		// update tdp
		tdp.process(currentSample);
		// log some stats
		final int tSize = currentSample.topology().nLeaves();
		outMan.write(treeName, 
				"treeSize", tSize, 
				"acceptPr", acceptPr, 
				//        "maskSparsity", currentSparsity,
				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
				"trans2tranv", trans2tranv, 
				"LogLikelihood", previousLogLLEstimate);
	}
}
