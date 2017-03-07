package pty.smc;


import ev.to.Ninja;
import fig.basic.Pair;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import pty.RootedTree;
import pty.smc.models.LikelihoodModelCalculator;
import ma.SequenceType;
import nuts.math.Sampling;
import nuts.util.Arbre;
import static nuts.util.CollUtils.*;

public class SimpleCoalescentState
{ 
	// public PartialCoalescentState clone;
	//	private  LikelihoodModelCalculator topLikelihoodModel=null;	
	//	private  double topBrLenU=-1.0, topBrLenL=-1.0;
	private  List<String> newNodeNames=list(); 
	private  double remainH; 
	private List<RootedTree> forest = list();
	private List<LikelihoodModelCalculator> modelCalculators = list();
	//	private Map<UnorderedPair<Integer, Integer>,Double> pairsToMerge=CollUtils.map();
	//	private List<Integer> pairToMerge= list();


	public void setRemainH(double newRemainH)
	{
		this.remainH=newRemainH; 
	}

	public double  getRemainH()
	{
		return remainH;
	}

	public void setNewNodesNames(List<String> newNodeNames)
	{
		this.newNodeNames=newNodeNames;
	}

	public String getNodesName(int i)
	{
		return newNodeNames.get(i); 
	}

	public String getNodesName()
	{
		//		System.out.println(newNodeNames);
		if(newNodeNames.size()==0) return ("internal_"+(nLeaves()-forest.size()));
		//		System.out.println("nRoots: "+nRoots()+"; size: "+newNodeNames.size()+" --->"+(newNodeNames.get(newNodeNames.size()-nRoots()+1)));
		return newNodeNames.get(newNodeNames.size()-nRoots()+1); 
	}

	//	public double getTopBrLenU()
	//	{
	//		return topBrLenU;
	//	}
	//
	//	public double getTopBrLenL()
	//	{
	//		return topBrLenL;
	//	}
	//
	//	public LikelihoodModelCalculator getTopLikelihoodModel()
	//	{
	//		return topLikelihoodModel;
	//	}

	public List<LikelihoodModelCalculator> getModelCalculators()
	{
		return modelCalculators;
	}


	//	public boolean noMessageFromTop()
	//	{
	//		return topLikelihoodModel==null || topBrLenU==-1.0|| topBrLenL==-1.0; 
	//	}

	public int indexOf(Taxon tax)
	{
		for(int i=0;i<nRoots();i++)
			if(forest.get(i).topology().getContents().equals(tax))return i;				
		return -1;
	}

	//	public static List<Pair<SimpleCoalescentState,Double>>  logWeightsOfFixedParticle(ParticleKernel<SimpleCoalescentState> kernel,RootedTree rt)
	//	{
	//		System.out.println("true tree: "+rt);
	//		Map<Pair<Taxon,Taxon>,Pair<Double,Double>> map=SimpleCoalescentState.restoreSequence(rt);
	//		return logWeights(kernel,map);
	//	}

	//	public static List<Pair<SimpleCoalescentState,Double>> logWeights(ParticleKernel<SimpleCoalescentState> kernel, 
	//			Map<Pair<Taxon,Taxon>,Pair<Double,Double>> sortMap)
	//	{     
	//		List<Pair<SimpleCoalescentState,Double>> result=list();
	//		SimpleCoalescentState current= kernel.getInitial();
	//		//initial
	//		Iterator<Pair<Taxon,Taxon>> arbreIterator=sortMap.keySet().iterator();
	//		while(arbreIterator.hasNext())
	//		{
	//			Pair<Taxon,Taxon> currentTaxaPair=arbreIterator.next();
	//			Pair<Double,Double> brLenPair=sortMap.get(currentTaxaPair); 
	//			Taxon first=currentTaxaPair.getFirst(),second=currentTaxaPair.getSecond();
	//			System.out.println("Pair to merge: "+ first+"	"+second);
	//			System.out.println("heights: "+ brLenPair.getFirst()+"	"+brLenPair.getSecond());
	//
	//			SimpleCoalescentState coalesceResult = current.coalesce(current.indexOf(first), 
	//					current.indexOf(second), brLenPair.getFirst(), brLenPair.getSecond());
	//			double logWeight = coalesceResult.logLikelihood()-current.logLikelihood();
	//			result.add(Pair.makePair(coalesceResult, logWeight));			
	//		} 
	//		return result;
	//	}


	//	//TODO:move this to a better place	
	//	public static List<Pair<SimpleCoalescentState,Double>> logWeights(ParticleKernel<SimpleCoalescentState> kernel, Map<Arbre<Taxon>,Double> sortHeightMap)
	//	{     
	//		List<Pair<SimpleCoalescentState,Double>> result=list();
	//		SimpleCoalescentState current= kernel.getInitial();
	//		//initial
	//		Iterator<Arbre<Taxon>> arbreIterator=sortHeightMap.keySet().iterator();
	//		double currentHeight=0;
	//		while(arbreIterator.hasNext())
	//		{
	//			Arbre<Taxon> currentArbre=arbreIterator.next();
	//			Arbre<Taxon> first=currentArbre.getChildren().get(0),second=currentArbre.getChildren().get(1);
	//			//			double newHeight=sortHeightMap.get(currentArbre);
	//			//			double delta=newHeight-currentHeight;
	//			//			System.out.println("heights: "+ newHeight+" - "+currentHeight+"="+delta);
	//			//			currentHeight=newHeight;			
	//			//			SimpleCoalescentState coalesceResult = current.coalesce(current.indexOf(first.getContents()), 
	//			//					current.indexOf(second.getContents()), delta);
	//			System.out.println("Pair to merge: "+ first.getContents()+"	"+second.getContents());
	//			System.out.println("heights: "+ sortHeightMap.get(first)+"	"+sortHeightMap.get(second));
	//			SimpleCoalescentState coalesceResult = current.coalesce(current.indexOf(first.getContents()), 
	//					current.indexOf(second.getContents()), sortHeightMap.get(first), sortHeightMap.get(second));
	//
	//			double logWeight = coalesceResult.logLikelihood()-current.logLikelihood();
	//			result.add(Pair.makePair(coalesceResult, logWeight));			
	//		} 
	//		return result;
	//	}



	//TODO: move this function to a better place	
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


	//TODO: move this function to a better place; use other ways to replace sortByValue?
	public static Pair<List<Pair<SimpleCoalescentState,Double>>,List<String>> restoreSequence(ParticleKernel<SimpleCoalescentState> kernel,RootedTree rt, double maxHeight, boolean uniformProp)
	{
		//double rtHeight=RootedTree.Util.height(rt);
		//		System.out.println("fixed rt is : "+rt);
		List<String> newNodeNames=list();
		List<Pair<SimpleCoalescentState,Double>> result=list();
		SimpleCoalescentState current= kernel.getInitial();
		List<Arbre<Taxon>> childrenList=rt.topology().nodes();
		//		List<UnorderedPair<Taxon,Taxon>> pairTaxaList =list();
		Map<Taxon,Double> branchLengths=rt.branchLengths();
		Map<Arbre<Taxon>,Double> heightMap=map();
		for(int i=0;i<childrenList.size();i++)
		{
			if(!childrenList.get(i).isLeaf()) 
				heightMap.put(childrenList.get(i), height(branchLengths, childrenList.get(i))); 		
		}
		//		for(Arbre<Taxon> arbre:heightMap.keySet())
		//		{
		//			//		System.out.println(arbre.deepToString());
		//			System.out.println("Height: "+heightMap.get(arbre));
		//		}
		//   sort the height in an increasing order		
		Map<Arbre<Taxon>,Double> sortHeightMap=sortByValue(heightMap);

		//		for(Arbre<Taxon> currentArbre:sortHeightMap.keySet())
		//		{
		//			System.out.println(currentArbre.getContents()+": "+sortHeightMap.get(currentArbre));
		//		}

		Set<Arbre<Taxon>> arbreSet=sortHeightMap.keySet();		
		Iterator<Arbre<Taxon>> arbreIterator=arbreSet.iterator();

		while(arbreIterator.hasNext())
		{			
			Arbre<Taxon> currentArbre=arbreIterator.next();
			newNodeNames.add(currentArbre.toString());
		}

		arbreIterator=arbreSet.iterator();
		//		Map<Pair<Taxon,Taxon>,Pair<Double,Double>> result=map();
		while(arbreIterator.hasNext())
		{			
			Arbre<Taxon> currentArbre=arbreIterator.next();
			Taxon first=currentArbre.getChildren().get(0).getContents(),second=currentArbre.getChildren().get(1).getContents();
			//						result.put(Pair.makePair(first, second), Pair.makePair(branchLengths.get(first),branchLengths.get(second)));
			//			current.printRoots();
			//			System.out.println("Pair to merge: "+ first+"	"+second);
			//			System.out.println("Indices: ("+current.indexOf(first)+", "+ current.indexOf(second)+")");             			
			//			System.out.println("heights: "+ branchLengths.get(first)+"	"+branchLengths.get(second));
			SimpleCoalescentState coalesceResult = current.coalesce(currentArbre.toString(),current.indexOf(first), 
					current.indexOf(second), branchLengths.get(first), branchLengths.get(second));
			coalesceResult.setNewNodesNames(newNodeNames);
			double newRemainH=maxHeight-sortHeightMap.get(currentArbre);
			coalesceResult.setRemainH(newRemainH);
			double logWeight = coalesceResult.logLikelihood()-current.logLikelihood();
			if(uniformProp)
			{
				double rate =1.0;
				int forestSize=current.nRoots();
				double lambda = rate*forestSize*(forestSize-1)*0.5;
				double delta=coalesceResult.height()-current.height();
				logWeight+=Sampling.exponentialLogDensity(1.0/lambda, delta)+Math.log(newRemainH);
			}
			current=coalesceResult;
			result.add(Pair.makePair(coalesceResult, logWeight));			
		}			


		//		System.out.println("restored tree: "+current.getForest().get(0));
		return Pair.makePair(result, newNodeNames);
	}


	/*
	public double gamma()
	{
		double result = 0.0;   		
		List<Integer> indexList=list();
		for(UnorderedPair<Integer, Integer> pair:pairsToMerge.keySet())
		{
		int first=pair.getFirst(),second=pair.getSecond();
		indexList.add(first);
		indexList.add(second);
		double currentBr=pairsToMerge.get(pair);
		SimpleCoalescentState newModelCal=coalesce(first, second, currentBr*0.5, currentBr*0.5); 
		result+=newModelCal.logLikelihood();		
		}

		for (int i=0; i<modelCalculators.size();i++)
		    if(!indexList.contains(i))
			result += modelCalculators.get(i).logLikelihood();
		return result;
	}
	 */

	/*	public UnrootedTree urtFromNj()
	{		 		
		Counter<UnorderedPair<Taxon, Taxon>> distances=new Counter<UnorderedPair<Taxon, Taxon>>(); 		
		for (int i = 0; i < nRoots(); i++)
			for (int j = i+1; j< nRoots(); j++)
			{
				Taxon t1 = forest.get(i).topology().getContents(), t2 = forest.get(j).topology().getContents();
				UnorderedPair<Taxon,Taxon> key = new UnorderedPair<Taxon,Taxon>(t1, t2);
				Pair<Double,Double> ss = DiscreteModelCalculator.k2pDistanceSuffStat((DiscreteModelCalculator) modelCalculators.get(i), (DiscreteModelCalculator) modelCalculators.get(j));
				distances.setCount(key, Ninja.k2pDistance(ss.getFirst(), ss.getSecond()));
			}
		return Ninja.inferTreeFromDistances(distances);
	}
	 */

	/*
	public UnrootedTree urtFromNj()
	{		 				
		List<Taxon> taxonList=list();
		List<Integer> indexList=list();

		for(UnorderedPair<Integer, Integer> pair:pairsToMerge.keySet())
		{
			indexList.add(pair.getFirst());
			taxonList.add(forest.get(pair.getFirst()).topology().getContents());		
			indexList.add(pair.getSecond());
			taxonList.add(forest.get(pair.getSecond()).topology().getContents());					
		}
		Counter<UnorderedPair<Taxon, Taxon>> distances=new Counter<UnorderedPair<Taxon, Taxon>>();

		for (int i = 0; i < indexList.size(); i++)
			for (int j = i+1; j< indexList.size(); j++)
			{
				UnorderedPair<Taxon,Taxon> key = new UnorderedPair<Taxon,Taxon>(taxonList.get(i), taxonList.get(j));
				Pair<Double,Double> ss = DiscreteModelCalculator.k2pDistanceSuffStat((DiscreteModelCalculator) modelCalculators.get(indexList.get(i)), (DiscreteModelCalculator) modelCalculators.get(indexList.get(j)));
				distances.setCount(key, Ninja.k2pDistance(ss.getFirst(), ss.getSecond()));
			}
		UnrootedTree urt=Ninja.inferTreeFromDistances(distances);
		Set<Taxon> vertexSet=urt.getTopology().vertexSet();		
		Iterator<Taxon> iter=vertexSet.iterator();
		List<Taxon> internalNodes=list();
		while(iter.hasNext())
		{
			Taxon currentTaxon=iter.next(); 			
			if(!taxonList.contains(currentTaxon)) // internal nodes
				internalNodes.add(currentTaxon); 
		}				
	   int nInternalNodes=internalNodes.size(); // either 1 or 2.
	   if(nInternalNodes==1) internalNodes.get(0); 

	   return urt;

	}	
	 */

	public double brFromNj(int first, int second)
	{		 				
		Pair<Double,Double> ss = DiscreteModelCalculator.k2pDistanceSuffStat((DiscreteModelCalculator) modelCalculators.get(first), (DiscreteModelCalculator) modelCalculators.get(second));
		return Ninja.k2pDistance(ss.getFirst(), ss.getSecond());
	}


	public static SimpleCoalescentState constructSimpleCoalescentState(List<RootedTree> forest, List<LikelihoodModelCalculator> modelCalculators)
	{
		SimpleCoalescentState result = new SimpleCoalescentState();		
		result.forest=forest;
		result.modelCalculators=modelCalculators;
		return result;
	}

	public int nRoots() { return forest.size(); }




	public static double height(RootedTree rt)
	{        		
		int nLeaves=rt.topology().nLeaves();
		double totalSum =0;
		double []pathLength=new double[nLeaves];
		Iterator<Arbre<Taxon>> leafIer = rt.topology().leaves().iterator();
		for(int i=0;i<nLeaves;i++)
		{
			double sum = 0;
			Arbre<Taxon> cur = leafIer.next();
			while (!cur.isRoot())
			{				
				//	System.out.println("is root? "+cur.isRoot()+" Branch length: "+rt.branchLengths().get(cur.getContents()));
				sum += rt.branchLengths().get(cur.getContents());
				cur = cur.getParent();          
			}
			pathLength[i]=sum;
			totalSum+=sum;
		}
		totalSum=totalSum/nLeaves;            
		return totalSum;
	}

	public static double sumBr(RootedTree rt)
	{        	  
		Map<Taxon,Double>  bl = rt.branchLengths();
		List<Arbre<Taxon>> nodes=rt.topology().nodes();
		double sum=0; 
		for(int i=0; i<nodes.size();i++)
		{
			//System.out.println(nodes.get(i).getContents());
			//System.out.println(bl.get(nodes.get(i).getContents()));
			if(!nodes.get(i).isRoot())
				sum+=bl.get(nodes.get(i).getContents()); 
		}	                
		return sum/rt.topology().nLeaves();
		//	return sum;
	}


	public double height()
	{
		double max = Double.NEGATIVE_INFINITY, current;
		for (RootedTree tree : forest)
			//  if ((current = RootedTree.Util.height(tree)) > max)
			if ((current = height(tree)) > max)
				//	 if ((current = sumBr(tree)) > max)
				max = current;

		return max;
	}




	public double subTreeHeightIncrease(int i)
	{

		//  System.out.println((height()-RootedTree.Util.height(forest.get(i))) +"	"+(height()-SimpleCoalescentState.height(forest.get(i))));
		// 	return height()-RootedTree.Util.height(forest.get(i));    	
		return height()-SimpleCoalescentState.height(forest.get(i));
		//		return height()-sumBr(forest.get(i));
	}


	public String toNewick()
	{
		return RootedTree.Util.toNewick(forest.get(0));
	}

	/* public double logLikelihood()
  {
    double result = 0.0;    
    for (LikelihoodModelCalculator calc : modelCalculators)     
      result += calc.logLikelihood();    
  //  System.out.println("number of trees: "+nRoots());
    //if(result!=clone.logLikelihood()) 
   // 	System.out.println(result+"-"+clone.logLikelihood()+"="+(result-clone.logLikelihood())); 
    return result;
  }
	 */



	public double logLikelihood()
	{
		double result = 0.0;   		
		for (int i=0; i<modelCalculators.size();i++)		
			result += modelCalculators.get(i).logLikelihood();		 
		return result;
	}



	public double weightedLikelihood()
	{		
		double result = 0.0;   		
		for (int i=0; i<modelCalculators.size();i++)
		{
			result += modelCalculators.get(i).extendLogLikelihood(subTreeHeightIncrease(i));
		}
		return result;
	}



	public double extendLogLikelihood2()
	{
		double result = 0.0;   		
		for (int i=0; i<modelCalculators.size();i++)
		{					
			//subTreeHeightIncrease(i)/height(forest.get(i));
			double tmp= height(forest.get(i))==0? modelCalculators.get(i).extendLogLikelihood(subTreeHeightIncrease(i)): modelCalculators.get(i).logLikelihood()*(1+height(forest.get(i))*0.002);
			result += tmp;			
		}		
		return result;
	}




	public double extendLogLikelihood()
	{
		double result = 0.0;   
		//System.out.println("Forest size: "+modelCalculators.size());		
		for (int i=0; i<modelCalculators.size();i++)
		{			
			//		System.out.println("Height increase: "+subTreeHeightIncrease(i)+"; Extended likelihood: "+modelCalculators.get(i).extendLogLikelihood(subTreeHeightIncrease(i))+"	; Likelihood: "+modelCalculators.get(i).logLikelihood()+"	");
			double tmp=modelCalculators.get(i).extendLogLikelihood(subTreeHeightIncrease(i));
			result += tmp;			
			//	System.out.println("Likelihood: "+tmp+" ; Height increase: "+subTreeHeightIncrease(i)+" ; Forest height: "+height()+" ;  subtree height:"+height(forest.get(i)));
		}		
		//  System.out.println("number of trees: "+nRoots());
		//if(result!=clone.logLikelihood()) 
		// 	System.out.println(result+"-"+clone.logLikelihood()+"="+(result-clone.logLikelihood())); 
		return result;
	}


	public double penalizedLogLikelihood(double alpha)
	{
		return logLikelihood()*Math.exp(alpha*(nRoots()-1));
	}



	public int nNontrivialSubtrees()
	{
		int count =0;
		for (RootedTree tree : forest)
			if (!tree.topology().isLeaf())
				count +=1;
		return count;
	}


	public static SimpleCoalescentState getInit(Map<Taxon,LikelihoodModelCalculator> likelihoodModels)
	{		
		SimpleCoalescentState result = new SimpleCoalescentState();
		for (Taxon t : likelihoodModels.keySet())
		{
			Arbre<Taxon> topo = Arbre.arbre(t);
			Map<Taxon,Double> bls = map();
			RootedTree rt = RootedTree.Util.create(topo, bls);
			result.forest.add(rt);
			result.modelCalculators.add(likelihoodModels.get(t));
		}    		
		return result;
	}


	//	public static SimpleCoalescentState getInit(Map<Taxon,LikelihoodModelCalculator> likelihoodModels,
	//			LikelihoodModelCalculator topLikelihoodModel, double topBrLenU, double topBrLenL)	
	//	{		
	//		SimpleCoalescentState result = new SimpleCoalescentState();
	//		for (Taxon t : likelihoodModels.keySet())
	//		{
	//			Arbre<Taxon> topo = Arbre.arbre(t);
	//			Map<Taxon,Double> bls = map();
	//			RootedTree rt = RootedTree.Util.create(topo, bls);
	//			result.forest.add(rt);
	//			result.modelCalculators.add(likelihoodModels.get(t));
	//		}    
	//		result.topLikelihoodModel=topLikelihoodModel;
	//		result.topBrLenU=topBrLenU;
	//		result.topBrLenL=topBrLenL;
	//		return result;
	//	}



	public static SimpleCoalescentState getInit(File alignment, SequenceType sequenceType)
	{
		Map<Taxon,DiscreteModelCalculator> inits = DiscreteModelCalculator.getInit(alignment, SequenceType.RNA);
		SimpleCoalescentState result = new SimpleCoalescentState();
		for (Taxon t : inits.keySet())
		{
			Arbre<Taxon> topo = Arbre.arbre(t);
			Map<Taxon,Double> bls = map();
			RootedTree rt = RootedTree.Util.create(topo, bls);
			result.forest.add(rt);
			result.modelCalculators.add(inits.get(t));
		}    
		//Dataset dataset = DatasetUtils.fromAlignment(alignment, SequenceType.RNA);
		//CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
		//   result.clone = PartialCoalescentState.initState(dataset, ctmc);
		//  System.out.println("init loglike "+result.logLikelihood()+"-"+result.clone.logLikelihood()+"="+(result.logLikelihood()-result.clone.logLikelihood())); 

		return result;
	}


	public List<RootedTree> getForest()
	{
		return forest;
	}


	public int nLeaves()
	{
		int count=0;
		for(int i=0;i<forest.size();i++)
		{
			count+=forest.get(i).topology().nLeaves();
		}
		return count; 
	}





	public SimpleCoalescentState coalesce(int left, int right, double delta)
	{
		// TODO: see RootedTree.Util.coalesce,  RootedTree.Util.height, and LikelihoodModelCalculator			
		double height=height();
		//		System.out.println("("+left+","+right+") "+delta);
		double leftHeight  = RootedTree.Util.height(forest.get(left));
		double rightHeight = RootedTree.Util.height(forest.get(right));
		double leftBranch=delta+height-leftHeight;
		double rightBranch=delta+height-rightHeight;
		return coalesce(left,right,  leftBranch, rightBranch);
	}

	public SimpleCoalescentState coalesce(String newNodeName, int left, int right, double delta)
	{
		// TODO: see RootedTree.Util.coalesce,  RootedTree.Util.height, and LikelihoodModelCalculator			
		double height=height();
		//				System.out.println(newNodeName+" ("+left+","+right+") "+delta);
		double leftHeight  = RootedTree.Util.height(forest.get(left));
		double rightHeight = RootedTree.Util.height(forest.get(right));
		double leftBranch=delta+height-leftHeight;
		double rightBranch=delta+height-rightHeight;
		SimpleCoalescentState newState=coalesce(newNodeName, left,right,  leftBranch, rightBranch);
		//		newState.topLikelihoodModel=this.topLikelihoodModel;
		newState.newNodeNames=this.newNodeNames;
		//		newState.topBrLenL=this.topBrLenL;
		//		newState.topBrLenU=this.topBrLenU;	
		return newState;
	}



	public SimpleCoalescentState coalesce(int left, int right, double leftBranchLength, double rightBranchlength)
	{
		return coalesce("internal_"+(nLeaves()-forest.size()), left,right, leftBranchLength, rightBranchlength);
	}

	public SimpleCoalescentState coalesce(String newNodeName, int left, int right, double leftBranchLength, double rightBranchlength)
	{

		//		System.out.println("Taxon name "+newNodeName+ " left "+forest.get(left).topology().getContents()+" right "+forest.get(right).topology().getContents()+" left branch: "+ leftBranchLength+" right branch: "+rightBranchlength);
		RootedTree newRoot = RootedTree.Util.coalesce(new Taxon(newNodeName), forest.get(left),forest.get(right), leftBranchLength, rightBranchlength);
		//RootedTree newRoot = RootedTree.Util.coalesce(new Taxon("internal_"+(nLeaves()-forest.size())), forest.get(left),forest.get(right), leftBranchLength, rightBranchlength);
		SimpleCoalescentState result=new SimpleCoalescentState();

		LikelihoodModelCalculator leftCal = modelCalculators.get(left);
		LikelihoodModelCalculator rightCal = modelCalculators.get(right);
		LikelihoodModelCalculator newCal=leftCal.combine(leftCal, rightCal, leftBranchLength, rightBranchlength, false);   

		for(int i=0;i<forest.size();i++)
		{
			if(i!=left && i!=right)
			{
				result.forest.add(forest.get(i)); 
				result.modelCalculators.add(modelCalculators.get(i));
			}
		}

		result.forest.add(newRoot);   
		result.modelCalculators.add(newCal);
		//  result.clone=clone.coalesce(left, right, delta, 0, 0);
		//  if(result.nRoots()!=this.nRoots()-1) throw new RuntimeException();
		return result;
	}


	public void printRoots()
	{
		System.out.print("nRoots are: ");
		for(int i=0;i<nRoots();i++)
		{
			System.out.print("("+i+": "+forest.get(i).topology().getContents()+"),	");
		}
		System.out.println();
	}

	@Override
	public String toString()
	{
		StringBuilder result= new StringBuilder();
		for (RootedTree rt : forest)
			result.append(rt + "\n");
		return result.toString();
	}


}
