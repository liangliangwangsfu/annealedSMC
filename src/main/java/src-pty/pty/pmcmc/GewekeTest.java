package pty.pmcmc;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import dr.math.distributions.NormalDistribution;
import pepper.Encodings;
import pty.RandomRootedTrees;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.RootedTree.RootingInfo;
import pty.io.Dataset;
import pty.io.Dataset.DatasetUtils;
import pty.smc.LazyParticleFilter;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.PriorPriorKernel;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.models.CTMC;
import ev.ex.DataGenerator;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import goblin.Taxon;
import ma.MSAPoset;
import ma.SequenceType;
import nuts.io.IO;
import nuts.math.Sampling;
import nuts.util.Counter;
import nuts.util.Indexer;

public class GewekeTest implements Runnable{
	@Option public int randseed=1;
	@Option public int nTaxa=5; 
	@Option public double treeRate=10.0;
	@Option public int M1=1000;
	@Option public int M2=1000;
	@Option public int len=50;
	@Option public int nParticles=10000;
	@Option public boolean isClock=true;
	@Option public double[] sigLevels=new double[]{0.05,.01,0.005,0.001};
	public static void main(String[] args)
	{
		IO.run(args,new GewekeTest(),"nc", NCPriorPriorKernel.class);
	}	

	@Override
	public void run() {		
		severalTestsPriorAndPriorPosterior(1,sigLevels);
		//		severalTestsPriorAndPrior(100,0.05);		
//				checkPrior(new Random(randseed),6,treeRate);
	}

	public int numberOfUrtTopo(int nTaxa)
	//The total number of distinct labelled topologies of an unrooted tree of n leaves is (2(n-1)-3)!!	
	{
		int result=1, i=(2*(nTaxa-1)-3);		
		while(i>1){
			result*=i;
			i-=2; 
		}
		LogInfo.logsForce("Taxa number is: "+nTaxa+"; The total number of distinct labelled topologies of an unrooted tree: "+result);
		return result;
	}	

	public int numberOfRtTopo(int nTaxa)
	//The total number of distinct labelled topologies of an unrooted tree of n leaves is (2(n-1)-3)!!	
	{
		int result=1, i=(2*nTaxa-3);		
		while(i>1){
			result*=i;
			i-=2; 
		}
		LogInfo.logsForce("Taxa number is: "+nTaxa+"; The total number of distinct labelled topologies of a rooted tree: "+result);
		return result;
	}	
	
	public Counter<UnrootedTree> getUniqueUnrootedTreeTopo(Random rand,int nTaxa,double treeRate)
	{		
		int nTrees=1000; 
		int noOfUrtTopo=numberOfUrtTopo(nTaxa);
		Counter<UnrootedTree> urtTopoSimu=simuUnrootedTreeTopo(rand,nTrees,nTaxa,treeRate);
		while(urtTopoSimu.size()!=noOfUrtTopo) 
		{
			nTrees*=2;
			urtTopoSimu=simuUnrootedTreeTopo(rand,nTrees,nTaxa,treeRate);
		}
		return urtTopoSimu;
	}

	public Counter<UnrootedTree> simuUnrootedTreeTopo(Random rand,int nTrees,int nTaxa,double treeRate)
	{		
		Counter<UnrootedTree> urtCounter = new Counter<UnrootedTree>();		
		for (int i = 0; i < nTrees; i++)
		{				  
			//			RootedTree  tree=TreeGenerators.sampleExpNonclock(rand, nTaxa, treeRate);
			RootedTree  tree=RandomRootedTrees.sampleCoalescent(rand,nTaxa,treeRate);
			UnrootedTree urt=UnrootedTree.fromRooted(tree);	
			Set<Set<Taxon>>  cladeSets=UnrootedTree.fromRooted(tree).unRootedClades();
			boolean newadd=true; 
			for(UnrootedTree key:urtCounter.keySet())
			{
				if(cladeSets.equals(key.unRootedClades()))
				{
					urtCounter.incrementCount(key, 1.0);
					newadd=false; 
					break; 
				}
			}
			if(newadd) 
			{
				UnrootedTree urtNoBranch=new UnrootedTree(urt);
				for(UnorderedPair<Taxon,Taxon> edge:urt.edges()) urtNoBranch.changeBranchLength(edge, 1); 
				urtCounter.incrementCount(urtNoBranch, 1.0);
			}
		}		  
		urtCounter.normalize();
		return urtCounter;
	}




	public void checkPrior(Random rand,int nTaxa,double treeRate)
	{			  
		//		Counter<UnrootedTree> urtCounter=simuUnrootedTreeTopo(rand,nTrees,nTaxa,treeRate);
		Counter<UnrootedTree> urtCounter=getUniqueUnrootedTreeTopo(rand,nTaxa,treeRate);
		int treeIndex=0;
//		int fourCounts=0,fiveCounts=0;
		double unifVar=0,firstMom=0,secondMom=0; 
		LogInfo.logsForce("Number of taxa: "+nTaxa+"; Number of simulated trees: "+urtCounter.size());
		for(UnrootedTree key:urtCounter.keySet())
		{			
			LogInfo.logsForce("Tree index: "+treeIndex++);
			LogInfo.logsForce(key.toNewick());		
			double diameter=key.diameter();
			firstMom+=diameter;
			secondMom+=Math.pow(diameter,2);
//			if(diameter==4.0) fourCounts++;
//			if(diameter==5.0) fiveCounts++;			
			LogInfo.logsForce("Probability: "+urtCounter.getCount(key)+" diameter: "+diameter);	
		}	
		int uniqueCount=urtCounter.size();
		unifVar=secondMom/uniqueCount-Math.pow((firstMom/uniqueCount), 2);
		LogInfo.logsForce("Mean and Variance of uniform tree topology: "+(firstMom/uniqueCount)+", "+unifVar);
//		LogInfo.logsForce("Four counts: "+fourCounts+"; five counts: "+fiveCounts);
	}

	public void severalTestsPriorAndPriorPosterior(int totalTest,double[] sigLevel)
	{
	
		SummaryStatistics pValuesBranch=new SummaryStatistics(),pValuesDiameter=new SummaryStatistics();
		Random rand=new Random(randseed);
		int[] sigcount1=new int[sigLevel.length],sigcount2=new int[sigLevel.length];
		for(int i=0;i<totalTest;i++){
			int randseed=Math.abs(rand.nextInt());		
			LogInfo.logsForce("\n--------Run "+(i+1)+"--------\n");
			Pair<Double,Double> pvalues=testPriorAndPriorPosterior(randseed);
			double pvalue1=pvalues.getFirst(),pvalue2=pvalues.getSecond();
			pValuesBranch.addValue(pvalue1);
			pValuesDiameter.addValue(pvalue2);
			for(int j=0;j<sigLevel.length;j++)
			{
			if(pvalue1<sigLevel[j]) sigcount1[j]++;
			if(pvalue2<sigLevel[j]) sigcount2[j]++;
			}
		}		
		LogInfo.logsForce("\n--------- Prior VS Prior-Posterior Summary: --------------"); 
		LogInfo.logsForce("Random number seed: "+randseed);
		LogInfo.logsForce("Number of taxa: "+nTaxa);
		LogInfo.logsForce("Sequence length: "+len); 
		LogInfo.logsForce("Number of particles: "+nParticles);
		LogInfo.logsForce("Clock tree? "+isClock);
		LogInfo.logsForce("Number of iterations: M1 "+M1+", M2 "+M2);
		LogInfo.logsForce("Total test runs: "+totalTest);
		LogInfo.logsForce("Testing the total branch length: ");
		LogInfo.logsForce("Significance levels 		Number of failed tests");
		for(int j=0;j<sigLevel.length;j++) LogInfo.logsForce(sigLevel[j]+"		"+sigcount1[j]);
		LogInfo.logsForce("Testing the tree diameter: ");
		LogInfo.logsForce("Significance levels 		Number of failed tests");
		for(int j=0;j<sigLevel.length;j++) LogInfo.logsForce(sigLevel[j]+"		"+sigcount2[j]);
	}


	public void severalTestsPriorAndPrior(int totalTest,double sigLevel)
	{
		SummaryStatistics pValuesBranch=new SummaryStatistics(),pValuesDiameter=new SummaryStatistics();
		Random rand=new Random(randseed);
		int sigcount1=0,sigcount2=0;
		for(int i=0;i<totalTest;i++){
			int randseed=rand.nextInt();		
			LogInfo.logsForce("\n--------Run "+(i+1)+"--------\n");
			Pair<Double,Double> pvalues=testPriorAndPrior(randseed);
			double pvalue1=pvalues.getFirst(),pvalue2=pvalues.getSecond();
			pValuesBranch.addValue(pvalue1);
			pValuesDiameter.addValue(pvalue2);
			if(pvalue1<sigLevel) sigcount1++;
			if(pvalue2<sigLevel) sigcount2++;
		}		
		LogInfo.logsForce("\n--------- Prior VS Prior Summary: --------------"); 
		LogInfo.logsForce("Random number seed: "+randseed);
		LogInfo.logsForce("Number of taxa: "+nTaxa);
		LogInfo.logsForce("Sequence length: "+len);
		LogInfo.logsForce("Clock tree? "+isClock);
		LogInfo.logsForce("Number of iterations: "+M1);
		LogInfo.logsForce("For testing the total branch length, "+sigcount1+" out of "+totalTest+" is less than "+sigLevel+";");
		LogInfo.logsForce("For testing the tree diameter, "+sigcount2+" out of "+totalTest+" is less than "+sigLevel+".");		
	}


	public Pair<Double,Double> testPriorAndPrior(int randseed){
		Random rand= new Random(randseed);
		Pair<SummaryStatistics,SummaryStatistics>  sumPrior=simuFromPrior(rand,nTaxa,treeRate,M1,isClock);
		Pair<SummaryStatistics,SummaryStatistics>  sumPrior2=simuFromPrior(rand,nTaxa,treeRate,M1,isClock);  
		double priormean1=sumPrior.getFirst().getMean(),priorvariance1=sumPrior.getFirst().getVariance();		//prior
		double priormean2=sumPrior2.getFirst().getMean(),priorvariance2=sumPrior2.getFirst().getVariance();   //posterior
		double stat=(priormean1-priormean2)/Math.sqrt(priorvariance1/M1+priorvariance2/M2);        
		NormalDistribution d=new NormalDistribution(0,1); 
		double pvalue=2*d.cdf(-Math.abs(stat));
		LogInfo.logsForce("Number of taxa: "+nTaxa);
		LogInfo.logsForce("Length of sequences: "+len);
		LogInfo.logsForce("Number of particles: "+nParticles);
		LogInfo.logsForce("Clock tree? "+isClock);
		LogInfo.logsForce("Number of iterations: M1 "+M1+", M2 "+M2);
		LogInfo.logsForce("Random number seed: "+randseed);

		LogInfo.logsForce("Total branch length: ");
		LogInfo.logsForce("Prior 1 --- Mean: "+priormean1+"; Variance: "+priorvariance1);
		LogInfo.logsForce("Prior 2 --- Mean: "+priormean2+"; Variance: "+priorvariance2);
		LogInfo.logsForce("Statistic "+stat+"; p value: "+pvalue);
		//
		priormean1=sumPrior.getSecond().getMean();
		priorvariance1=sumPrior.getSecond().getVariance();		//prior
		priormean2=sumPrior2.getSecond().getMean();
		priorvariance2=sumPrior2.getSecond().getVariance(); //posterior
		stat=(priormean1-priormean2)/Math.sqrt(priorvariance1/M1+priorvariance2/M2);                 
		double pvalue2=2*d.cdf(-Math.abs(stat));
		LogInfo.logsForce("Tree diameter: ");
		LogInfo.logsForce("Prior --- Mean: "+priormean1+"; Variance: "+priorvariance1);
		LogInfo.logsForce("Posterior --- Mean: "+priormean2+"; Variance: "+priorvariance2);
		LogInfo.logsForce("Statistic "+stat+"; p value: "+pvalue2);
		return Pair.makePair(pvalue, pvalue2);
	}

	public Pair<Double,Double> testPriorAndPriorPosterior(int randseed) {
		Random rand= new Random(randseed);
		Pair<SummaryStatistics,SummaryStatistics>  sumPrior=simuFromPrior(rand,nTaxa,treeRate,M1,isClock);
		Pair<SummaryStatistics,SummaryStatistics>  sumPosterior=simuFromPosterior(rand,nTaxa,len,treeRate,M2,nParticles,isClock);
		double priormean1=sumPrior.getFirst().getMean(),priorvariance1=sumPrior.getFirst().getVariance();		//prior
		double posteriormean1=sumPosterior.getFirst().getMean(),posteriorvariance1=sumPosterior.getFirst().getVariance(); //posterior
		double stat=(priormean1-posteriormean1)/Math.sqrt(priorvariance1/M1+posteriorvariance1/M2);        
		NormalDistribution d=new NormalDistribution(0,1); 
		double pvalue=2*d.cdf(-Math.abs(stat));
		LogInfo.logsForce("Number of taxa: "+nTaxa);
		LogInfo.logsForce("Length of sequences: "+len);
		LogInfo.logsForce("Number of particles: "+nParticles);
		LogInfo.logsForce("Clock tree? "+isClock);
		LogInfo.logsForce("Number of iterations: M1 "+M1+", M2 "+M2);
		LogInfo.logsForce("Random number seed: "+randseed);
		
		LogInfo.logsForce("Total branch length: ");
		LogInfo.logsForce("Prior --- Mean: "+priormean1+"; Variance: "+priorvariance1);
		LogInfo.logsForce("Prior-Posterior --- Mean: "+posteriormean1+"; Variance: "+posteriorvariance1);
		LogInfo.logsForce("Statistic "+stat+"; p value: "+pvalue);
		//
		priormean1=sumPrior.getSecond().getMean();
		priorvariance1=sumPrior.getSecond().getVariance();		   //prior
		posteriormean1=sumPosterior.getSecond().getMean();
		posteriorvariance1=sumPosterior.getSecond().getVariance();  //posterior
		stat=(priormean1-posteriormean1)/Math.sqrt(priorvariance1/M1+posteriorvariance1/M2);                 
		double pvalue2=2*d.cdf(-Math.abs(stat));
		LogInfo.logsForce("Tree diameter: ");
		LogInfo.logsForce("Prior --- Mean: "+priormean1+"; Variance: "+priorvariance1);
		LogInfo.logsForce("Posterior --- Mean: "+posteriormean1+"; Variance: "+posteriorvariance1);
		LogInfo.logsForce("Statistic "+stat+"; p value: "+pvalue2);
		return Pair.makePair(pvalue, pvalue2);
	}


	public  Pair<SummaryStatistics,SummaryStatistics> simuUniformUnrootedtree(Random rand,int nTaxa,double treeRate,int M1)
	{
		Counter<UnrootedTree> urtCounter=getUniqueUnrootedTreeTopo(rand,nTaxa,treeRate); 
		List<UnrootedTree> urtList= new ArrayList<UnrootedTree>(urtCounter.keySet());	    		
		SummaryStatistics brSum=new SummaryStatistics(),diameterSum=new SummaryStatistics();
		for(int i=0;i<M1;i++)
		{
			int randInd=rand.nextInt(urtCounter.size());
			//			System.out.println(randInd);
			UnrootedTree randTree=urtList.get(randInd);               
			UnrootedTree urt=new UnrootedTree(randTree);
			for(UnorderedPair<Taxon,Taxon> edge:randTree.edges()) urt.changeBranchLength(edge, Sampling.sampleExponential(rand,  1.0/treeRate));
			brSum.addValue(urt.totalBranchLength());           // add total tree length
			diameterSum.addValue(urt.diameter()); 			
		}
		return Pair.makePair(brSum, diameterSum);
	}


	public RootedTree sampleExpNonclock(Random rand,int nTaxa,double treeRate)
	{
		Counter<UnrootedTree> urtCounter=getUniqueUnrootedTreeTopo(rand,nTaxa,treeRate); 
		List<UnrootedTree> urtList= new ArrayList<UnrootedTree>(urtCounter.keySet());	    		
		int randInd=rand.nextInt(urtCounter.size());
		UnrootedTree randTree=urtList.get(randInd);               
		UnrootedTree urt=new UnrootedTree(randTree);
		for(UnorderedPair<Taxon,Taxon> edge:randTree.edges()) urt.changeBranchLength(edge, Sampling.sampleExponential(rand,  1.0/treeRate));
		UnorderedPair<Taxon,Taxon> randEdge=urt.randomEdge(rand); 
		RootingInfo rooting = new RootingInfo(randEdge.getFirst(), randEdge.getSecond(),new Taxon("internal_" + (nTaxa-2)),rand.nextDouble());
		return urt.reRoot(rooting); 		
	}


	public Pair<SummaryStatistics,SummaryStatistics> simuFromPrior(Random rand,int nTaxa,double treeRate,int M1,boolean isClock)
	{
		if(!isClock)return simuUniformUnrootedtree(rand,nTaxa,treeRate,M1);
		SummaryStatistics brSum=new SummaryStatistics(),diameterSum=new SummaryStatistics();
		for(int i=0;i<M1;i++)
		{
			RootedTree tree=null;
			tree=RandomRootedTrees.sampleCoalescent(rand, nTaxa, treeRate); // generate one random tree			
			UnrootedTree urt=UnrootedTree.fromRooted(tree);
			brSum.addValue(urt.totalBranchLength());           // add total tree length
			diameterSum.addValue(urt.diameter()); 			
		}
		return Pair.makePair(brSum, diameterSum);
	}

	public Pair<SummaryStatistics,SummaryStatistics>  simuFromPosterior(Random rand,int nTaxa,int len,double treeRate,int M2,int nParticles,boolean isPMCMC4clock)
	{		
		SummaryStatistics brSum=new SummaryStatistics(),diameterSum=new SummaryStatistics();
		RootedTree currentTree=null;
		if(isPMCMC4clock)
			currentTree=RandomRootedTrees.sampleCoalescent(rand,nTaxa,treeRate); // generate one random tree
		else  
			currentTree=sampleExpNonclock(rand,nTaxa,treeRate); 
		Indexer<Character> indexer = Encodings.dnaEncodings().nonGapCharactersIndexer();
		ParticleFilterOptions options = new ParticleFilterOptions();
		options.nParticles = nParticles;  
		options.nThreads = 1;   //TODO: solve the problems of using multiple threads in pmmh. 
		options.resampleLastRound = false;
		options.parallelizeFinalParticleProcessing = true;
		options.finalMaxNUniqueParticles = Integer.MAX_VALUE;
		options.maxNUniqueParticles = Integer.MAX_VALUE;
		options.rand = rand;
		options.verbose = false;
		SequenceType sequenceType=SequenceType.DNA;
		//		Random mainRand=new Random(123);
		//		double[] subsRates=Dirichlet.sample(rand, new double[]{10, 10, 10,10,10,10});
		// six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
		double[] subsRates=new double[]{0.26,0.18,0.17,0.15,0.11,0.13};  // six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
		//		double[] statFreqs=Dirichlet.sample(rand, new double[]{10,10,10,10});
		// stationary state frequencies. pi_A, pi_C, pi_G, pi_T				
		double[] statFreqs=new double[]{0.3,0.2,0.2,0.3};  // stationary state frequencies. pi_A, pi_C, pi_G, pi_T
		//		double alpha=Sampling.nextDouble(rand,0.1,0.9);       // shape parameter in the Gamma distribution
		double alpha=0.5;
		double pInv=0;        // the proportion of invariant sites
		double a_alpha=1.5,a_statFreqs=100,a_subsRates=100,a_pInv=0.2;
		int nCategories=4;
		int dataRepeatN=nCategories; 
		if(pInv>0)dataRepeatN=nCategories+1;		
		//		CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs,subsRates,4,len,alpha,nCategories,pInv);
		TreeDistancesProcessor tdp = new TreeDistancesProcessor();
		TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
		boolean useTopologyProcessor=false,useNJinfo=false;		
		int i=0;												
		while(i<M2)
		{
			i++;
			if(i%1000==1)System.out.print(i+"-");		 
			UnrootedTree urt=currentTree.getUnrooted();
			brSum.addValue(urt.totalBranchLength()); // add total tree length
			diameterSum.addValue(urt.diameter());
			// generate data using forward simulation
			MSAPoset msa=DataGenerator.generateDataFromGTRIGamma(rand,currentTree,indexer,true,len,pInv,alpha,statFreqs,subsRates);  // generateDNAdata=true			
			// PMCMC using the new dataset
			Dataset dataset=DatasetUtils.fromAlignment(msa,sequenceType,dataRepeatN);			
			StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();
			CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs,subsRates,4,dataset.nSites(),alpha,nCategories,pInv);
			PartialCoalescentState init = PartialCoalescentState.initFastState(dataset,ctmc,isPMCMC4clock); 
			LazyParticleKernel kernel = (isPMCMC4clock?new PriorPriorKernel(init):new NCPriorPriorKernel(init));			
			LazyParticleFilter<PartialCoalescentState> pf = new LazyParticleFilter<PartialCoalescentState>(kernel, options);					
			double zHat=	pf.sample(pro);
			double marginalLoglike=zHat+init.logLikelihood();	

			LazyParticleKernel kernel0 = (isPMCMC4clock?new PriorPriorKernel(init):new NCPriorPriorKernel(init));			
			LazyParticleFilter<PartialCoalescentState> pf0 = new LazyParticleFilter<PartialCoalescentState>(kernel0, options);					
			double zHat0=pf0.sample(pro);
			double marginalLoglike0=zHat0+init.logLikelihood();			
			double logRatio = marginalLoglike - marginalLoglike0;
			double acceptPr = Math.min(1, Math.exp(logRatio));
			final boolean accept = Sampling.sampleBern(acceptPr, rand);
			if (accept)
			{
				PartialCoalescentState sampled=pro.sample(rand);  //TODO: there is a problem here when the nThread >1				
				currentTree=sampled.getFullCoalescentState();
				//				System.out.println(currentTree.toString());
				// set ll!
				//				previousLogLLEstimate=marginalLoglike;
			}
			//			pmmh.setFixParameters(false);			
			//			pmmh.setSaveTreesFromPMCMC(false);
			//			pmmh.next(rand);    // one MCMC move
			// get the tree: pmmh.getRootedTree() 			
			//			currentTree=pmmh.getRootedTree();						
		}				    
		System.out.println();
		if(useTopologyProcessor){
			Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter(); 			    			    
			LogInfo.logsForce("\n Number of unique unrooted trees: "+urtCounter.keySet().size()); 			    
			for(UnrootedTree urt:urtCounter.keySet())
			{
				LogInfo.logsForce(urt);
				LogInfo.logsForce(urtCounter.getCount(urt));			    	
			}			
		}
		return Pair.makePair(brSum, diameterSum);
	}



}
