package pty.pmcmc;
import java.io.File;
import java.util.*;
import pty.RootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.smc.LazyParticleFilter;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.PriorPriorKernel;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.models.CTMC;
import ev.ex.NJPState2;
import ev.ex.NJStateNonclockKernel;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import fig.basic.ListUtils;
import fig.exec.Execution;
import fig.prob.Dirichlet;
import gep.util.OutputManager;
import goblin.Taxon;
import ma.MSAPoset;
import ma.MSAPoset.Column;
import nuts.io.IO;
import nuts.math.Sampling;

public class PMMH4GTRIGamma
{
	private final Dataset dataset;
	private ParticleFilterOptions options=null;
	private final TreeDistancesProcessor tdp;
	private boolean useTopologyProcessor=false;
	private final TreeTopologyProcessor trTopo; 
	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
	private RootedTree currentSample = null;
	private double[] subsRates;         // six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
	private double[] statFreqs;         // stationary state frequencies. pi_A, pi_C, pi_G, pi_T
	private double alpha=0.5;           // shape parameter in the Gamma distribution
	private double pInv=0.2;            // the proportion of invariant sites 
	private double a_alpha=1.3;         // tuning parameter for alpha. 
	private double a_pInv=0.2;          // tuning parameter for pInv. 
	private double a_statFreqs=300;     // tuning parameter for statFreqs; 
	private double a_subsRates=200;     // tuning parameter for subsRates;	
	public static OutputManager outMan = new OutputManager();
	private int iter=0; 
	private int nCategories=4;
	private int treeCount=0;
	File output=new File(Execution.getFile("results")); 
	private String nameOfAllTrees="allTrees.trees";
	private boolean saveTreesFromPMCMC=false;
	private boolean processTree=false;
	private boolean isPMMH4Clock=true;
	private boolean useNJinfo=false;  
	private boolean fixParameters=false;

	public PMMH4GTRIGamma(Dataset dataset0,ParticleFilterOptions options,TreeDistancesProcessor tdp,boolean useTopologyProcessor,
			TreeTopologyProcessor trTopo, double[] subsRates, double[] statFreqs, double alpha, 
			double pInv,double a_alpha,double a_pInv,double a_statFreqs,double a_subsRates,int nCategories,
			boolean processTree,boolean isPMMH4Clock,boolean useNJinfo)
	{
		this.dataset=dataset0;
		this.options = options;
		this.tdp = tdp;
		this.useTopologyProcessor=useTopologyProcessor;
		if(useTopologyProcessor)
			this.trTopo=trTopo; 
		else this.trTopo=null; 
		this.subsRates=subsRates;
		this.statFreqs=statFreqs; 
		this.alpha=alpha; 
		this.pInv=pInv;
		this.a_alpha=a_alpha;
		this.a_pInv=a_pInv;
		this.a_statFreqs=a_statFreqs;
		this.a_subsRates=a_subsRates;
		this.nCategories=nCategories;
		this.processTree=processTree;
		this.isPMMH4Clock=isPMMH4Clock;
		this.useNJinfo=useNJinfo; 
	}

	public void setFixParameters(boolean fixParameters)
	{
		this.fixParameters=fixParameters;
	}

	public boolean getFixParameters()
	{
		return this.fixParameters;
	}

	public void setSaveTreesFromPMCMC(boolean saveTreesFromPMCMC)
	{
		this.saveTreesFromPMCMC=saveTreesFromPMCMC; 
	}

	public void setNameOfAllTrees(String nameOfAllTrees)
	{
		this.nameOfAllTrees=nameOfAllTrees; 
	}

	public String getNameOfAllTrees()
	{
		return this.nameOfAllTrees;
	}

	public boolean getSaveTreesFromPMCMC()
	{
		return this.saveTreesFromPMCMC;
	}

	public void setProcessTree(boolean processTree)
	{
		this.processTree=processTree; 
	}

	public double[] getSubsRates()
	{
		return subsRates;
	}

	public double[] getStateFreqs()
	{
		return statFreqs; 
	}

	public double getAlpha()
	{
		return alpha; 
	}

	public double getpInv()
	{
		return pInv; 
	}

	public RootedTree getRootedTree()
	{
		return currentSample;
	}

	public double[] proposeAlpha(Random rand, double low, double high){
		double [] result=new double[2]; 
		double scale=0,proposedAlpha=Double.MAX_VALUE; 		
		while(proposedAlpha<low || proposedAlpha>high){
			scale=Sampling.nextDouble(rand, 1.0/a_alpha, a_alpha);
			proposedAlpha=scale*alpha;		
		}		
		result[0]=scale; 
		result[1]=proposedAlpha;
		return result; 
	}

	public double proposePInv(Random rand, double low, double high){		
		double proposedPInv=Double.MAX_VALUE; 		
		while(proposedPInv<low || proposedPInv>high){
			proposedPInv=Sampling.nextDouble(rand, Math.max(0, pInv-a_pInv), Math.min(1, pInv+a_pInv));
		}				
		return proposedPInv; 
	}

	public static double[] proposeFromDirichlet(Random rand, double a, double[] rates){		
		double[] alphas = new double[rates.length];
		for (int i = 0; i < rates.length; i++)alphas[i]=a*rates[i];
		double[] result = Dirichlet.sample(rand, alphas); 
		return result;
	}

	public double logProposal(double scale, double[] rates){
		double[] alphas = new double[rates.length];
		for (int i = 0; i < rates.length; i++)alphas[i]=scale*rates[i];	      	        	    		
		return Dirichlet.logProb(alphas, ListUtils.sum(alphas), rates);
	}

	public void next(Random rand)
	{
		iter++;
		RootedTree previousSample = currentSample;	
		int type;
		if(pInv==0)type=rand.nextInt(3);
		else type=rand.nextInt(4);
		// pInv: Sliding window 
		double proposedPInv=pInv; 
		double proposedAlpha=alpha, scale=1; 
		double[] proposedsubsRates=subsRates;
		double[] proposedstatFreqs=statFreqs;
		if(!fixParameters)
		{
			if(type==0){
				// proposals:
				// alpha: multiplier				
				double[] propAlpha=proposeAlpha(rand, 0.05, 50);
				scale=propAlpha[0];
				proposedAlpha=propAlpha[1];
			}else///
				if(type==1){
					// rates of substitutions
					proposedsubsRates=proposeFromDirichlet(rand, a_subsRates, subsRates);    
				}
				else
					if(type==2)
						proposedstatFreqs=proposeFromDirichlet(rand, a_statFreqs, statFreqs);
					else
						proposedPInv=proposePInv(rand, 0, 1);
		}
		// sample from PF
		StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();		 
		double acceptPr = 0.0;
		try 
		{
			CTMC ctmc = new CTMC.GTRIGammaCTMC(proposedstatFreqs,proposedsubsRates,4,dataset.nSites(),proposedAlpha,nCategories,proposedPInv);
			PartialCoalescentState init = PartialCoalescentState.initFastState(dataset,ctmc,isPMMH4Clock); 
			LazyParticleKernel kernel = (isPMMH4Clock?new PriorPriorKernel(init):(useNJinfo?new NJStateNonclockKernel(NJPState2.init(init)):new NCPriorPriorKernel(init)));			
			LazyParticleFilter<PartialCoalescentState> pf = new LazyParticleFilter<PartialCoalescentState>(kernel, options);					
			double zHat=	pf.sample(pro);			
			// compute accept/reject
			double marginalLoglike=zHat+init.logLikelihood();
			double logRatio = marginalLoglike - previousLogLLEstimate;
			if(!fixParameters)
			{
				if(type==0)
					logRatio += Math.log(scale);
				else
					if(type==1)
						logRatio +=logProposal(a_subsRates, proposedsubsRates)-logProposal(a_subsRates, subsRates);
					else
						if(type==2)
							logRatio +=logProposal(a_statFreqs, proposedstatFreqs)-logProposal(a_statFreqs, statFreqs);
			}
			acceptPr = Math.min(1, Math.exp(logRatio));
			if (currentSample != null &&  Double.isInfinite(acceptPr))
				throw new RuntimeException();
			final boolean accept = Sampling.sampleBern(acceptPr, rand);
//			System.out.println(marginalLoglike +" - "+previousLogLLEstimate+" = "+logRatio+"." +accept);
			if (accept)
			{
				PartialCoalescentState sampled=pro.sample(rand);  //TODO: there is a problem here when the nThread >1				
				currentSample=sampled.getFullCoalescentState();				
				// set ll!
				previousLogLLEstimate=marginalLoglike;
				statFreqs=proposedstatFreqs;
				subsRates=proposedsubsRates; 
				alpha=proposedAlpha;
				pInv=proposedPInv;				
			}
		}
		catch (Exception e)
		{
			// total fail!
			throw new RuntimeException();
		}
		// update tdp
		if(processTree)tdp.process(currentSample);
		if(useTopologyProcessor) trTopo.process(currentSample);
		String stringOfTree=RootedTree.Util.toNewick(currentSample);
		treeCount++;
		if(saveTreesFromPMCMC)
		{
			String cmdStr="echo -n 'TREE tree_" + treeCount+"=' >>"+nameOfAllTrees;
			IO.call("bash -s",cmdStr,output);
			cmdStr="echo '" +stringOfTree + "' | sed 's/internal_[0-9]*_[0-9]*//g'"+ " | sed 's/leaf_//g'"+ " >> "+nameOfAllTrees;
			IO.call("bash -s",cmdStr,output);
		}
		// log some stats
		final int tSize = currentSample.topology().nLeaves();
		outMan.write("PMMH",
				"Iter", iter,
				"treeSize", tSize,
				"acceptPr", acceptPr, 
				//        "maskSparsity", currentSparsity,
				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
				"statFreqs1", statFreqs[0],
				"statFreqs2", statFreqs[1],
				"statFreqs3", statFreqs[2],
				"statFreqs4", statFreqs[3],				
				"subsRates1",subsRates[0],
				"subsRates2",subsRates[1],
				"subsRates3",subsRates[2],
				"subsRates4",subsRates[3],
				"subsRates5",subsRates[4],
				"subsRates6",subsRates[5],
				"alpha",alpha,
				"pInv",pInv,
				"LogLikelihood", previousLogLLEstimate);
	}

	public static  boolean[]  findInvariantSites(MSAPoset align)
	{
		boolean[] result=new boolean[align.columns().size()];
		int i=0;
		for(Column c:align.linearizedColumns())
		{
			Map<Taxon, Integer> map=c.getPoints();			
			Character cha=null;
			boolean invariant=true;
			for(Taxon tax:map.keySet())
			{
				Character current=align.charAt(c, tax);
				if(cha==null) 
					cha=current;
				if(!current.equals(cha))  
				{
					invariant=false;
					break;
				}
			}
			result[i++]=invariant;					
		}		
		return result; 
	}
}
