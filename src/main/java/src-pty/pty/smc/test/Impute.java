package pty.smc.test;

import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.Dataset.DatasetType;
import pty.learn.DiscreteBP;
import pty.smc.PartialCoalescentState;
import pty.smc.models.DiscreteModelCalculator;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.CognateId;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import pty.eval.SymmetricDiff;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.PriorPostKernel;
import pty.smc.ParticleFilter.MAPDecoder;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import pty.smc.models.NoLikelihoodModel;
import nuts.util.Arbre;
import nuts.util.Counter;
import nuts.util.Arbre.ArbreMap;

import nuts.math.GMFct;
import nuts.math.GMFctUtils;
import nuts.math.Graph;
import nuts.math.TabularGMFct;
import nuts.math.TreeSumProd;

public class Impute implements Runnable {
	@Option public String filename;
	@Option boolean usePriorPost = true;
	
	public static final int N = 5;
	public static final Random rand= new Random(1);
	public static int[][] observations;
	public static ArrayList<Taxon> leafNames;
	public static int nleaves; 
	public static int nsites;
	
	public static void main (String args []){
		Execution.run(args, new Impute(), "hgdp", HGDPDataset.class);
	}
	
	
//	private static PartialCoalescentState createInitialState (String filename) {
//		PhylipParser pp = new PhylipParser (filename);
//		nleaves =  pp.getNumberOfLeaves();
//		nsites = pp.getNumberofSites();
//		observations = pp.getObservations();
//		leafNames = pp.getLeafNames();
//		
//		double[][] rate = {{-0.5,0.5}, {0.5,-0.5}} ;
//		CTMC ctmc = new CTMC.SimpleCTMC (rate, nsites);
//		ArrayList<LikelihoodModelCalculator> dmcList = new ArrayList <LikelihoodModelCalculator> ();
//		
//		/*
//		leafNames = new ArrayList<Language> ();
//		
//		
//		observations = new int[nleaves][];
//		for (int i = 0 ; i < nleaves; i++) {
//			observations [i] = new int[nsites];
//			for (int j = 0; j < nsites; j++) {
//				observations [i][j] = i%2; 
//				if ( i % 3 == 0 && j % 2 == 0)
//					observations [i][j] = -1;
//			} 
//		}*/
//		
//		for (int i = 0 ; i < nleaves; i++) {
//			dmcList.add( DiscreteModelCalculator.observation (ctmc, observations[i]));
//		}
//		PartialCoalescentState pcs = PartialCoalescentState.initialState( dmcList, leafNames);
//		return pcs;
//	}
	
 public static PartialCoalescentState initCTMCGeneState () { //String filename) {
    Dataset data = DatasetType.HGDP.loadDataset();
    int nsites = data.nSites();
    double rateScalar = 1e-8;
    double[][] rate = {{-rateScalar,+rateScalar}, {+rateScalar,-rateScalar}} ;
    CTMC ctmc = new CTMC.SimpleCTMC (rate, nsites);
    return  PartialCoalescentState.initState(data,ctmc);
  }

  public void run()
  {
	
    PartialCoalescentState pcs = initCTMCGeneState(); //createInitialState (filename);
    ParticleKernel <PartialCoalescentState> ppk = usePriorPost? 
    		new PriorPostKernel (pcs) : new PriorPriorKernel (pcs);
    ParticleFilter.ParticleMapperProcessor<PartialCoalescentState,PartialCoalescentState> mapd = 
      ParticleFilter.ParticleMapperProcessor.saveCoalescentParticlesProcessor();
    ParticleFilter.bootstrapFilter ( ppk , mapd, N , rand );
    
    ArrayList < Counter<String>> marginalProbabilities = new ArrayList < Counter <String>> ();
    for (int i = 0 ; i < nsites; i++) {
    	marginalProbabilities.add(new Counter<String>());
    }
    
    int count = 0 ;
    for (PartialCoalescentState p: mapd.getCounter().keySet()) {
    	System.out.println ("Processing tree " + count++);
    	Arbre<Taxon> tree = p.getArbreAndBranchLengths().getFirst();
    	List<Arbre<Taxon>> leaves = tree.getLeaves();
    	
    	for (int i = 0 ; i < nsites; i++) {
    		GMFct<Taxon> gmf = DiscreteBP.posteriorMarginalTransitions(p, i);
    		Counter <String> probCounter = marginalProbabilities.get(i);
    		for (Arbre<Taxon> l:leaves) {
    			probCounter.incrementCount (l.toString(), gmf.get(l.getContents(),1));	
    		}
    	}    
    	    	
    }
    
    fillIn (marginalProbabilities);
    
/*    	
	System.out.println ("Final");
	edu.berkeley.nlp.util.PriorityQueue<PartialCoalescentState> pq  = mapd.getCounter().asPriorityQueue();
	Arbre<Language> tree = pq.peek().getArbreAndBranchLengths().getFirst();
	ArrayList<Arbre<Language>> leaves = tree.getLeaves();
	for (Arbre<Language> l:leaves) {
		System.out.print("\"" + l.toString() + "\"\t");
	}
	System.out.print("\n");

	for (int i = 0 ; i < nsites; i++) {
		Counter <String> probCounter = marginalProbabilities.get(i);
		for (Arbre<Language> l:leaves) {
			probCounter.setCount (l.toString(), probCounter.getCount(l.toString())/N);
			System.out.print(probCounter.getCount(l.toString()) + "\t");
		}
		System.out.print("\n");
	}
*/
	
	
   
 }
  
  public void fillIn ( ArrayList<Counter<String>> marginalProbabilities) {
	  
	  for (int i = 0 ; i < nleaves; i++ ) {		  
		  System.out.print( leafNames.get(i) +  "\t");
		  for (int j = 0 ; j < nsites; j++) {
			  Counter <String> probCounter = marginalProbabilities.get(j);
			  if (observations[i][j] == -1)  {
				  double prob = probCounter.getCount(leafNames.get(i).toString());
				  int filledin = rand.nextDouble() < prob?1:0;
				  System.out.print(filledin + "\t"); 
			  } else 
				  System.out.print(observations[i][j] + "\t");
		  }
		  
		  System.out.print("\n");
	  }
  }
	
}