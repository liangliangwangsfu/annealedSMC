package smcsampler;
import java.io.*;
import java.util.*;

import ev.poi.processors.TreeDistancesProcessor;
import ev.ex.PhyloSamplerMain;
import ev.ex.TreeGenerators;
import ev.ex.PhyloSamplerMain.SStempScheme;
import fig.basic.Option;
import goblin.Taxon;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.mcmc.UnrootedTreeState;
import pty.smc.test.PhyloHash.Phylo;
import pty.mcmc.PhyloSampler;
import pty.mcmc.ProposalDistribution;
import pty.mcmc.PhyloSampler.PhyloProcessor;
import pty.mcmc.PhyloSampler.PhyloProcessorAdaptor;
import pty.mcmc.PhyloSampler.PriorOptions;

import ma.MSAParser;
import ma.SequenceType;
import nuts.io.IO;
import nuts.maxent.SloppyMath;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;

public class phyloMCMC2 implements Runnable{
	
	  @Option public File alignmentInputFile = null;
	  @Option public SequenceType st = SequenceType.RNA;
	  @Option public Random initTreeRandom = new Random(10);
	  @Option public UnrootedTree initTree =  null;
	 
	  private PhyloSampler sampler = null;
	  public UnrootedTree consensus = null;
	  public TreeDistancesProcessor tdp = null;
	  
	  private double[] temperatureSchedule;  //in the decreasing order
	  private int nSamplesEachChain=20000; 
	  private double logZ=0;
	  @Option public boolean computeLogZUsingSteppingStone=false;
	  
	  @Override
	  public void run()
	  {
	    sampler = new PhyloSampler();
	    sampler.setOutputText(false);
	    List<Taxon> leaves = MSAParser.parseMSA(alignmentInputFile).taxa();
	    //UnrootedTree initTree = UnrootedTree.fromRooted(TreeGenerators.sampleCoalescent(initTreeRandom,leaves, false));
	    //RootedTree initTree = TreeGenerators.sampleCoalescent(initTreeRandom,leaves, false);
	    //UnrootedTreeState ncts = UnrootedTreeState.fromAlignment(UnrootedTree.fromRooted(initTree), alignmentInputFile, st);
	    UnrootedTreeState ncts = UnrootedTreeState.fromAlignment(initTree, alignmentInputFile, st);
	    sampler.init(ncts);
	    tdp = new TreeDistancesProcessor();
	    PhyloProcessor tdpa = 
	      new PhyloSampler.PhyloProcessorAdaptor(tdp);
	    sampler.setProcessors(Arrays.asList(tdpa));    
	    sampler.sampleManyTimes();
	    // consensus decode
	    consensus = tdp.getConsensus();
		//int nSampleUsed=(int) (nSamplesEachChain*0.8);
		//int nBurnin=nSamplesEachChain-nSampleUsed;    
/*	    if(computeLogZUsingSteppingStone){
	    	sampler.init(sampler.getCurrentState());
	    	for(int i=0;i<temperatureSchedule.length-1;i++)
	    {
	    	sampler=sampler.createHeatedVersion(temperatureSchedule[i]);    	
	    
	    	double [] logweights = new double [nSampleUsed];
	    	double tempDiff=temperatureSchedule[i+1]-temperatureSchedule[i];
	    	//new added, previous temperature setting is problemastic
	    	sampler.setTemperature(temperatureSchedule[i]);
	    	if(i == 0) {
	    		for(int j=0; j < nSampleUsed; j++) {
	    			Random r  = new Random();
	    			initTree = TreeGenerators.sampleExpNonclock(r,leaves, 10.0);
	    			//initTree = TreeGenerators.sampleCoalescent(r,leaves, false);
	    			ncts = UnrootedTreeState.fromAlignment(UnrootedTree.fromRooted(initTree), alignmentInputFile, st);
	    		    logweights[j] = tempDiff*ncts.getLogLikelihood();
	    		}
	    	}else {
	        	for(int j=0;j<nSamplesEachChain;j++)
	        	{
	        	sampler.sample();
	        	if(j>=nBurnin) {logweights[j-nBurnin]=tempDiff*sampler.logLikelihood();
	        	}
	        	
	        	}
	    	}
	    //	setLogZ(getLogZ() + (SloppyMath.logAdd(logweights)-Math.log(1.0*nSampleUsed)));  
	    	System.out.println(logZ);*/
//	    }
//	    }    
	  }
	  
	  public static void main(String[] args)
	  {
	    IO.run(args, new PhyloSamplerMain(),
	        "prop", ProposalDistribution.Util._defaultProposalDistributionOptions,
	        "phylo", PhyloSampler._defaultPhyloSamplerOptions,
	        "prior", PhyloSampler._defaultPriorOptions);
	  }
	   

}
