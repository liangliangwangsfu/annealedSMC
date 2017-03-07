package pty.smc;
import java.util.*;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.pmcmc.ParticleGibbs4GTRIGamma;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.models.LikelihoodModelCalculator;
import ev.poi.processors.TreeDistancesProcessor;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import goblin.Taxon;
import nuts.util.Arbre;
import static nuts.util.CollUtils.*;


public class PG2
{
	private final ParticleFilter<PartialCoalescentState> pf;
	private final NCPriorPriorKernel kernel;
	private final TreeDistancesProcessor tdp;
	private final PartialCoalescentState model;

	private UnrootedTree current = null;

	public PG2(ParticleFilter<PartialCoalescentState> pf,
			NCPriorPriorKernel kernel, TreeDistancesProcessor tdp,
			PartialCoalescentState model)
	{
		this.pf = pf;
		this.kernel = kernel;
		this.tdp = tdp;
		this.model = model;
	}

	public void next(Random rand)
	{
		PartialCoalescentState.alwaysComputeTopMessage = true;
		StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();
		// if not initalized, just do an unconditional PF sampling round
		if (current == null)
		{
			pf.N=1000;
			System.out.println("initialized: pf.N is : "+pf.N);

			pf.sample(kernel, pro);
			// sample from the samples of the pf
			PartialCoalescentState sampled = pro.sample(rand);
			current = UnrootedTree.fromRooted(sampled.getFullCoalescentState());

			//     current = UnrootedTree.fromRooted(TreeGenerators.sampleExpNonclock(rand,  pf.N, 10));	    
			//    	current =  Ninja.inferTreeFromDistances(NJPState.refreshPairwiseDistance(new Counter<UnorderedPair<Taxon, Taxon>>(),model));
		}
		else
		{
			// pick a pair of taxa at random      

			//System.out.println(current); 
			RootedTree subtr1=null, subtr2=null;
			UnorderedPair<Taxon,Taxon>  selectEdge=null;
			//    	 boolean trySelectSubtr=true;
			//    	 while(trySelectSubtr){
			selectEdge = current.randomEdge(rand);
			List<Taxon> leaveTaxa=current.leaves();

			while(leaveTaxa.contains(selectEdge.getFirst()) || leaveTaxa.contains(selectEdge.getSecond()))
				selectEdge = current.randomEdge(rand);

			//System.out.println(selectEdge); 

			Pair<RootedTree,RootedTree> twoSubRootedtrees=Unrooted2RootedTree.divideOneUnrootedTree2twoRootedTrees(current, selectEdge);      

			if(twoSubRootedtrees.getFirst().nTaxa()<twoSubRootedtrees.getSecond().nTaxa())
			{
				subtr1=twoSubRootedtrees.getFirst();      
				subtr2=twoSubRootedtrees.getSecond();
			}else{
				subtr1=twoSubRootedtrees.getSecond();      
				subtr2=twoSubRootedtrees.getFirst();    	  
			}
			//      if(subtr1.nTaxa()<current.nTaxa()*0.5) trySelectSubtr=false;
			//   	 }
			double selectEdgeBl=current.branchLength(selectEdge);
			pf.N=subtr1.nTaxa()*10000;
			System.out.println(" pf.N is : "+pf.N);

			//      System.out.println("First subtree: "); 
			//     System.out.println(subtr1);

			//      System.out.println("Second subtree: ");
			//      System.out.println(twoSubRootedtrees.getSecond());

			List<Taxon> leavesList=subtr1.topology().leaveContents();           
			PartialCoalescentState  init=init(leavesList);
			List<Pair<PartialCoalescentState,Double>> restorePCS=ParticleGibbs4GTRIGamma.restoreSequence4NonClockTree(init,subtr1);
			List<PartialCoalescentState> path=list();
			double[] weights=new double[restorePCS.size()];
			for(int i=0;i<restorePCS.size();i++){
				path.add(restorePCS.get(i).getFirst());			 
				weights[i]=restorePCS.get(i).getSecond();
				//			System.out.println(weights[i]);			
			}			 		



			//System.out.println(path);
			/*      for(int i=0;i<path.size();i++)
      {
    	  System.out.println(i+"-th element in the path: ");
    	  System.out.println(path.get(i).toString());
      }
			 */

			// set init state
			kernel.setInitial(init);

			/*      // compute weights
      double [] weights  = new double[path.size()];
      for (int i = 0; i < weights.length; i++)
        weights[i] =  (NCPriorPriorKernel.useOptimal ? path.get(i).nonClockLogWeight() : path.get(i).logLikelihoodRatio() - Math.log(path.get(i).nNonTrivialRoots()));
			 */    
			// set the conditioning and its weights
			pf.setConditional(path, weights);

			// do the sampling
			pf.sample(kernel, pro);

			// sample from the samples of the pf
			PartialCoalescentState sampled = pro.sample(rand);
			//      System.out.println("sampled: ");             
			//     System.out.println(UnrootedTree.fromRooted(sampled.getFullCoalescentState()));

			// reconstitute full tree
			//current = dissector.reconstitute(sampled); //sampled.getFullCoalescentState();
			RootedTree pooledTree=pool2SubTrees(subtr1, subtr2, selectEdge,  selectEdgeBl);
			//      System.out.println("pooled tree: ");             
			//     System.out.println(pooledTree);

			current=pooledTree.getUnrooted(); 
		}
		// process current sample!
		tdp.process(current);
	}


	private PartialCoalescentState  init(List<Taxon> taxaToConsider)
	{
		if (model.isClock())
			throw new RuntimeException();
		//    if (topology.getChildren().size() == 2)
		//     halfTopBranch = originalTree.branchLength(topology.getChildren().get(0).getContents(), topology.getChildren().get(1).getContents()) / 2.0;
		List<Taxon> leavesNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		for (Arbre<CoalescentNode> currentArbre : model.roots)
		{
			if (!currentArbre.isLeaf()) throw new RuntimeException();
			Taxon t = currentArbre.getContents().nodeIdentifier;
			if (taxaToConsider.contains(t))
			{
				leavesNames.add(t);
				leaves.add(currentArbre.getContents().likelihoodModelCache);
			}
		}
		return PartialCoalescentState.initialState(
				leaves, 
				leavesNames,
				model.getObservations(),
				false);    
	}


	//static public RootedTree mergeTwoSubtrees(RootedTree mainTree, RootedTree subTree, Taxon mergePoint)

	public RootedTree pool2SubTrees(RootedTree subTree1, RootedTree subTree2, UnorderedPair<Taxon,Taxon>  connectEdge, double selectEdgeBl)
	{

		// pool topology
		Arbre<Taxon> arbre1=subTree1.topology().copy(),arbre2=subTree2.topology().copy();
		final List<Arbre<Taxon>>  twoChildren=new ArrayList<Arbre<Taxon>>();
		twoChildren.add(arbre1);
		twoChildren.add(arbre2);
		Arbre<Taxon> arbre=Arbre.arbre(Taxon.dummy, twoChildren);

		// pool branch lengths. 
		Map<Taxon,Double> branchLengths=map();// brLens1=subTree1.branchLengths(),brLens2=subTree2.branchLengths();
		branchLengths.putAll(subTree1.branchLengths());
		branchLengths.putAll(subTree2.branchLengths());
		double halfBl=selectEdgeBl*0.5;
		branchLengths.put(subTree2.topology().getContents(),halfBl);
		branchLengths.put(subTree1.topology().getContents(),halfBl);		 
		return new RootedTree.Util.RootedTreeImpl(arbre,branchLengths);		
	}

}
