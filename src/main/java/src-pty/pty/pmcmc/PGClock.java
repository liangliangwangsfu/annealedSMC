package pty.pmcmc;


import static nuts.util.CollUtils.list;
import static nuts.util.CollUtils.map;

import java.util.*;

import pty.RootedTree;
import pty.io.Dataset;
import pty.io.TreeEvaluator;

import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.PriorPriorKernel;
import pty.smc.models.CTMC;
//import punt.SimpleCoalescentState;

import fig.basic.Pair;
import gep.util.OutputManager;
import goblin.Taxon;
import ev.poi.processors.TreeDistancesProcessor;
import nuts.math.Sampling;
import nuts.util.Arbre;


public class PGClock 
{
//	private final Dataset dataset;
//	private final ParticleFilter<PartialCoalescentState> pf;
//	private ParticleKernel<PartialCoalescentState> kernel=null;
//	private final TreeDistancesProcessor tdp;
//	private double previousLogLLEstimate = Double.NEGATIVE_INFINITY;
//	private RootedTree currentSample = null;
//	private double trans2tranv;
//	public double a=1.5;
//	public static OutputManager outMan = new OutputManager();
//	public boolean stop=false;
//	private int tryCountBeforAccept=0;
//	public String out="PGClock";
//	private int iter=0;
//	private PartialCoalescentState conditionedTree=null; 
//
//
//	public PGClock(Dataset dataset0, ParticleFilter<PartialCoalescentState> pf,  TreeDistancesProcessor tdp, double trans2tranv0)
//	{
//		this.trans2tranv=trans2tranv0;
//		this.dataset=dataset0;
//		this.pf = pf;
//		this.tdp = tdp;
//	}
//
//	public void next(Random rand)
//	{
//		iter++;
//		RootedTree previousSample = currentSample;
//		// an MH step to sample the parameter kappa.
//		double scale=Sampling.nextDouble(rand, 1.0/a, a);
//		double proposedTrans2tranv = scale*trans2tranv;
//		double proposedLoglike = 0, currentLoglike=0;   
//		final double logRatio =proposedLoglike  - currentLoglike + 0.5*trans2tranv*(1-scale) + Math.log(scale);
//		double acceptPr = Math.min(1, Math.exp(logRatio));
//		final boolean accept=Sampling.sampleBern(acceptPr, rand);
//		if (accept) trans2tranv=proposedTrans2tranv;				
//		
//		CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), trans2tranv);  
//		PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, ctmc, true);  // clock	          					
//		kernel = new PriorPriorKernel(init);
//
//   	 // I need to find a path that leads to the reconstruction of a given clock tree 
//	//	or a non-clock tree.
//				
//		// conditionedTree should be a clock tree.
//		if(!conditionedTree.isFinalState() || !conditionedTree.isClock()) throw new RuntimeException();  
//		conditionedTree.get
//	    List<PartialCoalescentState> path =  dissector.currentSampledPath(model, rand);
//		
//		conditionedTree.getArbreAndBranchLengths(); 
//		
//		
//		
//	      PartialCoalescentState init = path.get(0);
//	      path = path.subList(1, path.size());
//	      
//	      // set init state
//	      kernel.setInitial(init);
//	      
//	      // compute weights
//	      double [] weights  = new double[path.size()];
//	      for (int i = 0; i < weights.length; i++)
//	        weights[i] =  (NCPriorPriorKernel.useOptimal ? path.get(i).nonClockLogWeight() : path.get(i).logLikelihoodRatio() - Math.log(path.get(i).nNonTrivialRoots()));
//	      
//	      // set the conditioning and its weights
//	      pf.setConditional(path, weights);
//	      
//	      // do the sampling
//	      pf.sample(kernel, pro);
//	      
//	      // sample from the samples of the pf
//	      PartialCoalescentState sampled = pro.sample(rand);
//
//		
//		
//		
//		
//		
//		
//		StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();		 
//		double acceptPr = 0.0;
////		double[] currentMetricResults=new double[TreeEvaluator.coreTreeMetrics.size()];
////		double[] metricVec=new double[TreeEvaluator.coreTreeMetrics.size()];		
//		try 
//		{
//			tryCountBeforAccept++;
////			if(tryCountBeforAccept>10) stop=true;
//
//			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), proposedTrans2tranv);  
//			PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, ctmc, true);  // clock	          					
//			kernel = new PriorPriorKernel(init);
//			//			kernel.nonclockTreeRate=nonclockTreeRate;			
////			if(tryCountBeforAccept>5 && pf.N<10000)
////			{
////				pf.N=pf.N+500;
////				tryCountBeforAccept=0;
////			}
////			pf.sample(kernel, processors);
//			pf.sample(kernel, pro);
//			// compute accept/reject
//			double marginalLoglike=pf.estimateNormalizer()+init.logLikelihood();
//			final double logRatio =marginalLoglike  - previousLogLLEstimate + 0.5*trans2tranv*(1-scale) + Math.log(scale);
////						System.out.println("log(scale) "+Math.log(scale)+" proposed MarginalLike "+pf.estimateNormalizer()+" Currentlog MarginalLike: "+previousLogLLEstimate+"; log Ratio: "+logRatio);
//			acceptPr = Math.min(1, Math.exp(logRatio));
////			System.out.println(pf.N+" "+acceptPr+": "+marginalLoglike  +" -"+ previousLogLLEstimate+(0.5*trans2tranv*(1-scale)+Math.log(scale)));
//			if (currentSample != null &&  Double.isInfinite(acceptPr))
//				throw new RuntimeException();
////			UnrootedTree  inferred = innerTdp.getConsensus();				
////			for (int i=0; i<TreeEvaluator.coreTreeMetrics.size(); i++)
////			{
////				currentMetricResults[i] = TreeEvaluator.coreTreeMetrics.get(i).score(inferred, goldut);
////				metricVec[i]=acceptPr*currentMetricResults[i]+(1-acceptPr)*metricResults[i];
////			}
//
//			final boolean accept = Sampling.sampleBern(acceptPr, rand);
//			if (accept)
//			{
////				tryCountBeforAccept=0;
//				// sample from sample
//				PartialCoalescentState sampled = pro.sample(rand);				
//				currentSample = sampled.getFullCoalescentState();
//				// set ll!
//				previousLogLLEstimate = marginalLoglike;
//				trans2tranv=proposedTrans2tranv;
////				for (int i=0; i<TreeEvaluator.coreTreeMetrics.size(); i++) metricResults[i] = currentMetricResults[i];
//			}
//		}
//		catch (Exception e)
//		{
//			// total fail!
//		}
//
////				System.out.println(currentSample.topology());
//		// update tdp
//		tdp.process(currentSample);
//		// log some stats
//		final int tSize = currentSample.topology().nLeaves();
//		outMan.write(out,
//				"Iter", iter,
//				"treeSize", tSize, 
//				"pfN",pf.N, 
//				"acceptPr", acceptPr, 
//				//        "maskSparsity", currentSparsity,
//				"rfDist", (previousSample == null ? 0 : new TreeEvaluator.RobinsonFouldsMetric().score(currentSample, previousSample)),
//				"trans2tranv", trans2tranv, 
//				"LogLikelihood", previousLogLLEstimate
////				TreeEvaluator.coreTreeMetrics.get(0), metricVec[0],
////				TreeEvaluator.coreTreeMetrics.get(1), metricVec[1],
////				TreeEvaluator.coreTreeMetrics.get(2), metricVec[2],
////				TreeEvaluator.coreTreeMetrics.get(3), metricVec[3],
////				TreeEvaluator.coreTreeMetrics.get(4), metricVec[4]
//				);
//	}
//	
//	
//	public static Pair<List<Pair<PartialCoalescentState,Double>>,List<String>> restoreSequence(ParticleKernel<PartialCoalescentState> kernel, 
//			PartialCoalescentState conditioned)  
//	{
//		if(!conditioned.isFinalState() || !conditioned.isClock()) throw new RuntimeException();
//		List<PartialCoalescentState> conditional= list();  
//		List<String> newNodeNames=list();
//		double [] conditionalUnnormWeights=new double[conditioned.getLeaves()-1]; 
//		PartialCoalescentState current= kernel.getInitial();
//		RootedTree rt=conditioned.getFullCoalescentState(); 		
//		
//		List<Arbre<Taxon>> childrenList=rt.topology().nodes(); 
//		Map<Taxon,Double> branchLengths=rt.branchLengths();
//		Map<Arbre<Taxon>,Double> heightMap=map();
//		for(int i=0;i<childrenList.size();i++)
//		{
//			if(!childrenList.get(i).isLeaf()) 
//				heightMap.put(childrenList.get(i), height(branchLengths, childrenList.get(i))); 		
//		}
//		Map<Arbre<Taxon>,Double> sortHeightMap=sortByValue(heightMap);
//		Set<Arbre<Taxon>> arbreSet=sortHeightMap.keySet();		
//		Iterator<Arbre<Taxon>> arbreIterator=arbreSet.iterator();
//		while(arbreIterator.hasNext())
//		{			
//			Arbre<Taxon> currentArbre=arbreIterator.next();
//			newNodeNames.add(currentArbre.toString());
//		}
//		arbreIterator=arbreSet.iterator();
//		while(arbreIterator.hasNext())
//		{						
//			Taxon taxon=current.getRoots().get(0).getContents().nodeIdentifier; 
//			
//			
//			Arbre<Taxon> currentArbre=arbreIterator.next();
//			Taxon first=currentArbre.getChildren().get(0).getContents(),second=currentArbre.getChildren().get(1).getContents();
//			
//			int firstIndex=indexOf(current,first), secondIndex=indexOf(current,second); 
//			
//			PartialCoalescentState coalesceResult = current.coalesce(currentArbre.toString(), current.indexOf(first), 
//					current.indexOf(second), branchLengths.get(first), branchLengths.get(second));
//		
//			
//			
//			
//			
//			
//			
//			coalesceResult.setNewNodesNames(newNodeNames);
//			double newRemainH=maxHeight-sortHeightMap.get(currentArbre);
//			coalesceResult.setRemainH(newRemainH);
//			double logWeight = coalesceResult.logLikelihood()-current.logLikelihood();
//			if(uniformProp)
//			{
//				double rate =1.0;
//				int forestSize=current.nRoots();
//				double lambda = rate*forestSize*(forestSize-1)*0.5;
//				double delta=coalesceResult.height()-current.height();
//				logWeight+=Sampling.exponentialLogDensity(1.0/lambda, delta)+Math.log(newRemainH);
//			}
//			current=coalesceResult;
//			result.add(Pair.makePair(coalesceResult, logWeight));			
//		}			
//
//
//		//		System.out.println("restored tree: "+current.getForest().get(0));
//		return Pair.makePair(result, newNodeNames);
//	}
//	
//	public static double height(Map<Taxon,Double> branchLengths, Arbre<Taxon> arbre)
//	{
//		double sum=0; 
//		do{
//			arbre=arbre.getChildren().get(0);
//			sum+=branchLengths.get(arbre.getContents());         			   
//		}while(!arbre.isLeaf());
//
//		return sum;
//	}
//
//	public static Map sortByValue(Map map) {
//		List list = new LinkedList(map.entrySet());
//		Collections.sort(list, new Comparator() {
//			public int compare(Object o1, Object o2) {
//				return ((Comparable) ((Map.Entry) (o1)).getValue())
//				.compareTo(((Map.Entry) (o2)).getValue());
//			}
//		});
//
//		Map result = new LinkedHashMap();
//		for (Iterator it = list.iterator(); it.hasNext();) {
//			Map.Entry entry = (Map.Entry)it.next();
//			result.put(entry.getKey(), entry.getValue());
//		}
//		return result;
//	} 
//
//	public static int indexOf(PartialCoalescentState currentState, Taxon tax)
//	{		
//		
//		for(int i=0;i<currentState.nRoots();i++)
//			if(currentState.getRoots().get(i).getContents().equals(tax)) return i;				
//		return -1;
//	}
//
//
}
