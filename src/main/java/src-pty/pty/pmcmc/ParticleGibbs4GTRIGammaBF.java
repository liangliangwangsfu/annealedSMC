package pty.pmcmc;

import static nuts.util.CollUtils.list;
import static nuts.util.CollUtils.map;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import nuts.io.IO;
import nuts.math.Sampling;
import nuts.util.Arbre;
import pty.RootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;
import pty.mcmc.UnrootedTreeState;
import pty.smc.BackForwardKernel0;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.PartialCoalescentState;
import pty.smc.PartialCoalescentState4BackForwardKernel;
import pty.smc.ParticleFilter;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.ParticleKernel;
import pty.smc.models.CTMC;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import fig.basic.ListUtils;
import fig.basic.Pair;
import fig.exec.Execution;
import fig.prob.Dirichlet;
import gep.util.OutputManager;
import goblin.Taxon;

public class ParticleGibbs4GTRIGammaBF {
	private final Dataset dataset;
	ParticleFilterOptions options=null;
	private final TreeDistancesProcessor tdp;
	private boolean useTopologyProcessor=false;	
	private final TreeTopologyProcessor trTopo;
	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
	private RootedTree currentSample = null;
	private double[] subsRates;  // six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
	private double[] statFreqs;  // stationary state frequencies. pi_A, pi_C, pi_G, pi_T
	private double alpha=0;       // shape parameter in the Gamma distribution
	private double pInv=0;        // the proportion of invariant sites 
	private double a_alpha=1.3;       // tuning parameter for alpha. 
	private double a_pInv=0.4;        // tuning parameter for pInv. 
	private double a_statFreqs=300;     // tuning parameter for statFreqs; 
	private double a_subsRates=200;     // tuning parameter for subsRates;	
	public static OutputManager outMan = new OutputManager();		
	private int iter=0; 
	private int nCategories=4;
	private int treeCount=0;	
	File output=new File(Execution.getFile("results")); 
	private String nameOfAllTrees="allTrees.trees";
	private boolean saveTreesFromPMCMC=false;
	private int sampleTreeEveryNIter=100; 
	private boolean processTree=false;
	private boolean isGS4Clock=true;
	
	public ParticleGibbs4GTRIGammaBF(Dataset dataset0,
			ParticleFilterOptions options, TreeDistancesProcessor tdp,
			boolean useTopologyProcessor,TreeTopologyProcessor trTopo,
			RootedTree initrt, double[] subsRates, double[] statFreqs, double alpha, double pInv,double a_alpha,double a_pInv,double a_statFreqs,double a_subsRates,int nCategories,
			boolean processTree,boolean isGS4Clock,int sampleTreeEveryNIter)
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
		this.currentSample=initrt;   
		this.a_alpha=a_alpha;
		this.a_pInv=a_pInv;
		this.a_statFreqs=a_statFreqs;
		this.a_subsRates=a_subsRates;
		this.nCategories=nCategories;
		this.processTree=processTree;
		this.isGS4Clock=isGS4Clock;
		this.sampleTreeEveryNIter=sampleTreeEveryNIter;
	}

	public ParticleGibbs4GTRIGammaBF(Dataset dataset0,
			ParticleFilterOptions options, TreeDistancesProcessor tdp,
			RootedTree initrt, double[] subsRates, double[] statFreqs,
			double alpha, double pInv, TreeTopologyProcessor trTopo)
	{
		this.dataset=dataset0;
		this.options = options;
		this.tdp = tdp;
		this.subsRates=subsRates;
		this.statFreqs=statFreqs; 
		this.alpha=alpha; 
		this.pInv=pInv; 
		this.currentSample=initrt;
		this.trTopo=trTopo;
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


	public double[] proposeFromDirichlet(Random rand, double a, double[] rates){		
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
		// pInv: Sliding window 
		double proposedPInv=0; 
		if(pInv>0)proposedPInv=proposePInv(rand, 0, 1);
		double acceptpInv=MHpInv(proposedPInv,rand);
		// proposals:
		// alpha: multiplier				    
		double[] propAlpha=proposeAlpha(rand, 0.05, 50);
		double scale=propAlpha[0];
		double proposedAlpha=propAlpha[1];
		double acceptPralpha=MHalpha(proposedAlpha,scale,rand); 
		// rates of substitutions
		double[] proposedsubsRates=proposeFromDirichlet(rand,a_subsRates,subsRates);   
		double acceptPrsubsRates=MHsubsRates(proposedsubsRates,a_subsRates,rand);
		double[] proposedstatFreqs=proposeFromDirichlet(rand,a_statFreqs,statFreqs);
		double acceptPrstatFreqs=MHstatFreqs(proposedstatFreqs,a_statFreqs,rand);
		// sample from PF
		StoreProcessor<PartialCoalescentState4BackForwardKernel> pro = new StoreProcessor<PartialCoalescentState4BackForwardKernel>();
		if((iter % sampleTreeEveryNIter) == 0)
		{
			CTMC currentctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4, dataset.nSites(), alpha, nCategories, pInv);
			PartialCoalescentState init0 = PartialCoalescentState
					.initFastState(dataset, currentctmc, true);
			PartialCoalescentState4BackForwardKernel init = new PartialCoalescentState4BackForwardKernel(
					init0, null, 0);

			ParticleKernel<PartialCoalescentState4BackForwardKernel> kernel = new BackForwardKernel0(
					init);

			ParticleFilter<PartialCoalescentState4BackForwardKernel> pf = new ParticleFilter<PartialCoalescentState4BackForwardKernel>();
			pf.nThreads = 1;
			pf.resampleLastRound = false;
			List<Pair<PartialCoalescentState4BackForwardKernel, Double>> restorePCS = restoreSequence(
					kernel, currentSample, isGS4Clock);
			List<PartialCoalescentState4BackForwardKernel> path = list();
			double[] weights=new double[restorePCS.size()];
			for(int i=0;i<restorePCS.size();i++){
				path.add(restorePCS.get(i).getFirst());			 
				weights[i]=restorePCS.get(i).getSecond();
			}
			// set the conditioning and its weights
			pf.setConditional(path, weights);

			// do the sampling			
			pf.sample(kernel, pro);
			PartialCoalescentState sampled = pro.sample(rand);
			currentSample=sampled.getFullCoalescentState();
			// update tdp
			if(processTree)tdp.process(currentSample);
			if(useTopologyProcessor) trTopo.process(currentSample);
			treeCount++;
			if(saveTreesFromPMCMC)
			{
			String stringOfTree=RootedTree.Util.toNewick(currentSample);
			String cmdStr="echo -n 'TREE gsTree_" + treeCount+"=' >>"+nameOfAllTrees;
			IO.call("bash -s",cmdStr,output);
			cmdStr="echo '" +stringOfTree + "' | sed 's/internal_[0-9]*_[0-9]*//g'"+ " | sed 's/leaf_//g'"+ " >> "+nameOfAllTrees;
			IO.call("bash -s",cmdStr,output);
			}
		}
		// log some stats
		final int tSize = currentSample.topology().nLeaves();
		outMan.write("PGS",
				"Iter", iter,
				"treeSize", tSize,
				"acceptPralpha", acceptPralpha,		
				"acceptpInv", acceptpInv,		
				"acceptPrstatFreqs", acceptPrstatFreqs,
				"acceptPrsubsRates", acceptPrsubsRates,
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

	private void updateLogLikelihood(){
		CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4, dataset.nSites(), alpha, nCategories, pInv);
		UnrootedTreeState previousncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);
		previousLogLLEstimate=previousncs.logLikelihood();
	}	

	private double MHalpha(double proposedAlpha,double scale,Random rand){
		CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4, dataset.nSites(), proposedAlpha, nCategories, pInv);
		UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);    
		double logratio=ncs.logLikelihood()-previousLogLLEstimate+Math.log(scale);
		double acceptPr= Math.min(1, Math.exp(logratio)); 
		final boolean accept = Sampling.sampleBern(acceptPr, rand);
		if (accept)
		{
			alpha=proposedAlpha; 
			updateLogLikelihood();
		}
		return acceptPr;
	}

	private double MHpInv(double proposedpInv,Random rand){
		CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4, dataset.nSites(), alpha, nCategories, proposedpInv);
		UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);    
		double logratio=ncs.logLikelihood()-previousLogLLEstimate;
		double acceptPr= Math.min(1, Math.exp(logratio)); 
		final boolean accept = Sampling.sampleBern(acceptPr, rand);
		if (accept)
		{
			pInv=proposedpInv; 
			updateLogLikelihood();
		}
		return acceptPr;
	}

	private double MHstatFreqs(double[] proposedstatFreqs,double a_statFreqs,Random rand){
		CTMC ctmc = new CTMC.GTRIGammaCTMC(proposedstatFreqs, subsRates, 4, dataset.nSites(), alpha, nCategories, pInv);
		UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);
		double logratio=ncs.logLikelihood()-previousLogLLEstimate+logProposal(a_statFreqs, proposedstatFreqs)-logProposal(a_statFreqs, statFreqs);
		double acceptPr= Math.min(1, Math.exp(logratio)); 
		final boolean accept = Sampling.sampleBern(acceptPr, rand);
		if (accept)
		{
			statFreqs=proposedstatFreqs; 
			updateLogLikelihood();
		}
		return acceptPr;
	}

	private double MHsubsRates(double[] proposedsubsRates,double a_subsRates,Random rand){
		CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, proposedsubsRates, 4, dataset.nSites(), alpha, nCategories, pInv);
		UnrootedTreeState ncs = UnrootedTreeState.initFastState(currentSample.getUnrooted(), dataset, ctmc);
		double logratio=ncs.logLikelihood()-previousLogLLEstimate+logProposal(a_subsRates, proposedsubsRates)-logProposal(a_subsRates, subsRates);     		
		double acceptPr= Math.min(1, Math.exp(logratio)); 
		final boolean accept = Sampling.sampleBern(acceptPr, rand);
		if (accept)
		{    
			subsRates=proposedsubsRates; 
			updateLogLikelihood();
		}
		return acceptPr;
	}

	public static double height(Map<Taxon,Double> branchLengths, Arbre<Taxon> arbre)
	{
		double sum=0; 
		do{
			arbre=arbre.getChildren().get(0);
			sum+=branchLengths.get(arbre.getContents());         			   
		}while(!arbre.isLeaf());
		return sum;
	}

	public static Map sortByValue(Map map) {
		List list = new LinkedList(map.entrySet());
		Collections.sort(list, new Comparator() {
			public int compare(Object o1, Object o2) {
				return ((Comparable) ((Map.Entry) (o1)).getValue())
						.compareTo(((Map.Entry) (o2)).getValue());
			}
		});
		Map result = new LinkedHashMap();
		for (Iterator it = list.iterator(); it.hasNext();) {
			Map.Entry entry = (Map.Entry)it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	} 

	public static List<Pair<PartialCoalescentState,Double>> restoreSequence4NonClockTree(PartialCoalescentState current,RootedTree rt)
	{
		List<Pair<PartialCoalescentState,Double>> result=list();			 		
		List<Arbre<Taxon>> childrenList=rt.topology().nodes();
		class ArbreComparator implements Comparator<Arbre<Taxon>> {
			@Override
			public int compare(Arbre<Taxon> arbre1, Arbre<Taxon> arbre2) {		    	
				return   arbre1.getContents().toString().compareTo(arbre2.getContents().toString());      
			}
		}
		Collections.sort(childrenList, new ArbreComparator());
		Map<Taxon,Double> branchLengths=rt.branchLengths();
		for(int i=0;i<childrenList.size();i++)
		{
			if(!childrenList.get(i).isLeaf())
			{
				Arbre<Taxon> currentArbre=childrenList.get(i);
				Taxon first=currentArbre.getChildren().get(0).getContents(),second=currentArbre.getChildren().get(1).getContents();
				PartialCoalescentState coalesceResult = current.coalesce(current.indexOf(first), current.indexOf(second), 0, branchLengths.get(first), branchLengths.get(second),currentArbre.getContents());							
				double logWeight = coalesceResult.logLikelihoodRatio()-Math.log(coalesceResult.nNonTrivialRoots());  //coalesceResult.logLikelihood()-current.logLikelihood();								
				current=coalesceResult;
				result.add(Pair.makePair(coalesceResult, logWeight));
			}
		}		
		return result;
	}

	public static List<Pair<PartialCoalescentState4BackForwardKernel, Double>> restoreSequence(
			ParticleKernel<PartialCoalescentState4BackForwardKernel> kernel,
			RootedTree rt, boolean isClock)
	{
		// if(!isClock)return
		// restoreSequence4NonClockTree(kernel.getInitial(),rt);
		List<String> newNodeNames=list();
		List<Pair<PartialCoalescentState4BackForwardKernel, Double>> result = list();
		PartialCoalescentState4BackForwardKernel current = kernel.getInitial();
		List<Arbre<Taxon>> childrenList=rt.topology().nodes();
		Map<Taxon,Double> branchLengths=rt.branchLengths();
		Map<Arbre<Taxon>,Double> heightMap=map();
		for(int i=0;i<childrenList.size();i++)
		{
			if(!childrenList.get(i).isLeaf()) 
				heightMap.put(childrenList.get(i), height(branchLengths, childrenList.get(i))); 		
		}
		//   sort the height in an increasing order		
		Map<Arbre<Taxon>,Double> sortHeightMap=sortByValue(heightMap);
		Set<Arbre<Taxon>> arbreSet=sortHeightMap.keySet();		
		Iterator<Arbre<Taxon>> arbreIterator=arbreSet.iterator();
		while(arbreIterator.hasNext())
		{			
			Arbre<Taxon> currentArbre=arbreIterator.next();
			newNodeNames.add(currentArbre.toString());
		}
		arbreIterator=arbreSet.iterator();
		double previousHeight=0;
		while(arbreIterator.hasNext())
		{			
			Arbre<Taxon> currentArbre=arbreIterator.next();
			Taxon first=currentArbre.getChildren().get(0).getContents(),second=currentArbre.getChildren().get(1).getContents();
			double currentHeight=sortHeightMap.get(currentArbre);
			double currentDelta=currentHeight-previousHeight;
			previousHeight=currentHeight;
			PartialCoalescentState4BackForwardKernel coalesceResult = (PartialCoalescentState4BackForwardKernel) current
					.coalesce(current.indexOf(first), current.indexOf(second),
							currentDelta, 0, 0, currentArbre.getContents());
			coalesceResult.setDeltaOld(currentDelta);
			coalesceResult.setParent(current);
			// PartialCoalescentState4BackForwardKernel;
			double logWeight = coalesceResult.logLikelihoodRatio();// coalesceResult.logLikelihood()-current.logLikelihood();
			current=coalesceResult;
			result.add(Pair.makePair(coalesceResult, logWeight));			
		}			
		return result;
	}
}



