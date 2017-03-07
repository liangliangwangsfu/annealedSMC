package pty.smc;
import java.util.*;

import pty.UnrootedTree;
import pty.smc.ParticleFilter.StoreProcessor;
import ev.poi.processors.TreeDistancesProcessor;
import goblin.Taxon;
import nuts.math.Sampling;

public class PG
{
	private final ParticleFilter<PartialCoalescentState> pf;
	private final NCPriorPriorKernel kernel;
	private final TreeDistancesProcessor tdp;
	private final PartialCoalescentState model;
	private UnrootedTree current = null;

	public PG(ParticleFilter<PartialCoalescentState> pf,
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
			pf.sample(kernel, pro);
			// sample from the samples of the pf
			PartialCoalescentState sampled = pro.sample(rand);
			current = UnrootedTree.fromRooted(sampled.getFullCoalescentState());

			//current = UnrootedTree.fromRooted(TreeGenerators.sampleExpNonclock(rand, 32, 10));	    
			//     current =  Ninja.inferTreeFromDistances(NJPState.refreshPairwiseDistance(new Counter<UnorderedPair<Taxon, Taxon>>(),model));
		}
		else
		{
			// pick a pair of taxa at random
			boolean success = false;
			CondSMCDissector dissector = null;
			//      System.out.println("current unrooted tree:");
			//      System.out.println(current.getUnrooted());
			subLoop:for (int i = 0; i < 100; i++)
			{
				List<Taxon> list = current.leaves();
				List<Integer> pair = Sampling.sampleWithoutReplacement(rand, list.size(), 2);
				Taxon t1 = list.get(pair.get(0)), t2 = list.get(pair.get(1));

				//        System.out.println("Two taxa: "+t1+", "+t2);
				// create the dissector
				dissector = new CondSMCDissector(current, t1, t2);
				if (dissector.sampleAnchor(rand))
				{
					success = true;
					break subLoop;
				}
			}
			if (!success)
				throw new RuntimeException();
			List<PartialCoalescentState> path = dissector.currentSampledPath(model, rand);
			PartialCoalescentState init = path.get(0);
			path = path.subList(1, path.size());
			//      System.out.println(path);
			/*      for(int i=0;i<path.size();i++)
      {
    	  System.out.println(i+"-th element in the path: ");
    	  System.out.println(path.get(i).toString());
      }
			 */    
			// set init state
			kernel.setInitial(init);

			// compute weights
			double [] weights  = new double[path.size()];
			for (int i = 0; i < weights.length; i++)
				weights[i] =  (NCPriorPriorKernel.useOptimal ? path.get(i).nonClockLogWeight() : path.get(i).logLikelihoodRatio() - Math.log(path.get(i).nNonTrivialRoots()));

			pf.N=weights.length*1000;
			System.out.println(" pf.N is : "+pf.N);

			// set the conditioning and its weights
			pf.setConditional(path, weights);

			// do the sampling
			pf.sample(kernel, pro);

			// sample from the samples of the pf
			PartialCoalescentState sampled = pro.sample(rand);       

			// reconstitute full tree
			current = dissector.reconstitute(sampled); //sampled.getFullCoalescentState();
			//      System.out.println(current);
		}
		// process current sample!
		tdp.process(current);
	}
}
