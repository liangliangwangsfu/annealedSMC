package smcsampler;

import static nuts.util.CollUtils.list;
import static nuts.util.CollUtils.union;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import ma.MSAParser;
import ma.SequenceType;
import nuts.io.CSV;
import nuts.io.IO;
import nuts.math.StatisticsMap.DescriptiveStatisticsMap;
import nuts.util.Counter;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.io.Dataset.DatasetUtils;
import pty.io.TreeEvaluator;
import pty.io.TreeEvaluator.TreeMetric;
import pty.mcmc.PhyloSampler;
//import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import pty.pmcmc.PhyloPFSchedule;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.PriorPriorKernel;
import pty.smc.models.CTMC;
import conifer.trees.StandardNonClockPriorDensity;
import ev.ex.PhyloSamplerMain;
import ev.ex.TreeGenerators;
import ev.poi.processors.TreeDistancesProcessor;
import ev.to.NJ;
import fig.basic.IOUtils;
import fig.basic.LogInfo;

import fig.basic.Option;
import fig.exec.Execution;
import fig.prob.Gamma;
import goblin.Taxon;

public class SMCSamplerExperiments implements Runnable
{
	@Option public boolean resampleRoot = false;
	@Option public InferenceMethod refMethod = InferenceMethod.MB;
	@Option public double nThousandIters = 10;
	@Option public ArrayList<InferenceMethod> methods = list(Arrays.asList(InferenceMethod.ANNEALING,InferenceMethod.MB));
	@Option public ArrayList<Double> iterScalings = list(Arrays.asList(1.0));
	@Option public int refIterScaling = 100;
	@Option public int repPerDataPt = 1;    
	@Option public double pmcmcSMCExpMix = 2.0/3.0;
	@Option public double treeRatePrior=1.0;
	public static ev.ex.DataGenerator.DataGeneratorMain generator = new ev.ex.DataGenerator.DataGeneratorMain();
	public static PhyloSamplerMain samplerMain = new PhyloSamplerMain();
	@Option 	public static Random mainRand  = new Random(3);
	@Option public boolean verbose = false;
	@Option public int nThreads = 1;
	@Option public int maxNUniqueParticles = Integer.MAX_VALUE;
	@Option public int finalMaxNUniqueParticles = Integer.MAX_VALUE;
	@Option public int nParticlesEachStep=1000; 
	@Option public double parameter_a=1.2;  
	@Option public double a_alpha=1.5; 
	@Option public double a_statFreqs=100;
	@Option public double a_subsRates=100;
	@Option public double a_pInv=0.2; 
	@Option public double pmmhPgsExpMix = 0.05;
	@Option public double burninPercent=0.2; // in PMCMC algorithm
	@Option public int sampleTreeEveryNIter=500;
	@Option public String RcommandDir="/Users/lwa68/workspace/alpha-liangliang/Rcode/";
	//@Option public String RcommandDir="/Users/l.wang/workspace/alpha-liangliang/Rcode/";
	@Option public boolean useDataGenerator=true; 	
	@Option public File dataFile=null;
	@Option public File  refTree=null;
	@Option public String dataDirName = "output";		
	@Option public SequenceType sequenceType=SequenceType.RNA; 
	@Option public boolean isPMCMC4clock=true;
	@Option 	public boolean  useNJinfo=false;	
	@Option 	public boolean  useTopologyProcessor=false; 
	@Option 	public boolean  saveTreesFromPMCMC=false;
	@Option 	public String nameOfAllTrees="allTrees.trees";
	@Option 	public double csmc_trans2tranv=2.0;
	@Option public double smcmcmcMix = 0.5;
	@Option public boolean betterStartVal=true;
	@Option public SMCSampler.ResamplingStrategy resamplingStrategy=SMCSampler.ResamplingStrategy.ESS;
	@Option
	public int nCSMC = 10;
	@Option
	public int nUCSMC = 10;
	@Option
	public boolean adaptiveTempDiff = false;
	@Option
	public int adaptiveType = 0;

	@Option
	public double alphaSMCSampler = 0.95;

	@Option
	public double essRatioThreshold = 0.7;
	@Option public int nNumericalIntegration = 100000;

	@Option
	public int nAnnealing = 5000;
	public  File data = null;
	private File output = null;
	private RootedTree goldrt;
	private PrintWriter logZout = null;
	private double marginalLogLike=0;

	public static void main(String[] args)
	{
		IO.run(args, new SMCSamplerExperiments(), 
				"pcs", PartialCoalescentState.class,
				"mb",   MrBayes.instance, 
				"gen", generator,
				"mcmc", samplerMain,
				"prop", ProposalDistribution.Util._defaultProposalDistributionOptions,
				"phylo", PhyloSampler._defaultPhyloSamplerOptions,
				"prior", PhyloSampler._defaultPriorOptions,
				"nj", NJ.class,
				"nc", NCPriorPriorKernel.class,
				"priorprior", PriorPriorKernel.class,
				"annealingKernel", AnnealingKernel.class
				);
	}





	protected void treeComparison()
	{
		//    data = new File( generator.output, "sim-0.msf");
		logZout = IOUtils.openOutEasy(new File(Execution.getFile("results"),
				"logZout.csv"));
		logZout.println(CSV.header("treeName", "Method", "NumericalLogZ", "logZ"));
		PrintWriter out = IOUtils.openOutEasy( new File(output, "results.csv"));
		out.println(CSV.header("Method", "Adaptive", "AdaptiveType", "T",
				"IterScale",
				"Repeat",
				"Metric", "Value", "TreeName", "Time"));
		List<File> files =null;	
		if(useDataGenerator)
			files = IO.ls(generator.output, "msf");
		else
			files = IO.ls(new File(Execution.getFile(dataDirName)), "msf");
		//    ReportProgress.progressBlock(files.size());

		boolean adaptiveTempDiff0=adaptiveTempDiff;
		int nMrBayesIter=0;		
		for (File f : files)
		{
			this.data = f;
			final String treeName = f.getName();
			LogInfo.track("Current tree:" + treeName);

			// evaluate the likelihood of the inferred tree
			Dataset dataset = DatasetUtils.fromAlignment(this.data, sequenceType);
			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());

			UnrootedTree goldut = 
					(generator.useGutellData ||!useDataGenerator)?
							(refTree==null?null:UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(refTree)))):
								UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(new File(f.getAbsolutePath().replaceAll("[.]msf$",".newick")))));

							File computedRefTrees = new File(Execution.getFile("computed-ref-trees"));
							//      ReportProgress.progressBlock(methods.size());
							for (int j = 0; j < repPerDataPt; j++)
							{
								for (int i = 0; i < methods.size(); i++) {
									InferenceMethod m = methods.get(i);

									int nRun = 1;
									if (m==InferenceMethod.ANNEALING && adaptiveTempDiff)
										nRun = 2;
									int nIter = Integer.MIN_VALUE;
									for (int l = 0; l < nRun; l++) {
										double iterScale = iterScalings.get(i);
										if(m == InferenceMethod.MB && nMrBayesIter>0) iterScale=nMrBayesIter;
										LogInfo.track("Current method:" + m + " with iterScale=" + iterScale + " (i.e. " + (iterScale * nThousandIters * 1000.0) + " iterations)");

										DescriptiveStatisticsMap<String> stats = new DescriptiveStatisticsMap<String>();
										//        ReportProgress.progressBlock(repPerDataPt);
										String treeNameCurrentRep= treeName;
										treeNameCurrentRep = treeNameCurrentRep + ".Rep"+ j;
										LogInfo.track("Repeat " + (j+1) + "/" + repPerDataPt);
										//          LogInfo.forceSilent = true;
										long  time=System.currentTimeMillis();
										TreeDistancesProcessor processor = m.doIt(this, iterScale, goldut, treeNameCurrentRep);									
										time=System.currentTimeMillis()-time;
										LogInfo.forceSilent = false;
										UnrootedTree inferred = processor.getConsensus();		
										System.out.println(inferred);
										IO.writeToDisk(new File(output, "consensus_"+treeName.replaceAll("[.]msf$",".newick")), inferred.toNewick());
										{
											// evaluate the likelihood of the inferred tree
											//											dataset = DatasetUtils.fromAlignment(this.data, sequenceType);
											//											ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
											UnrootedTreeState ncs = UnrootedTreeState.initFastState(inferred, dataset, ctmc);

											if (m == InferenceMethod.ANNEALING)
												out.println(CSV.body(m, adaptiveTempDiff, adaptiveType,
														this.nAnnealing, iterScale, j,
														"ConsensusLogLL", ncs.logLikelihood(),
														treeName, time));
											else
												out.println(CSV.body(m, "", "","", iterScale, j,
														"ConsensusLogLL", ncs.logLikelihood(),
														treeName, time));											
										}
										{
											// best log likelihood, when available
											double bestLogLL = processor.getBestLogLikelihood();
											if (m == InferenceMethod.ANNEALING)
												out.println(CSV.body(m, adaptiveTempDiff, adaptiveType,
														this.nAnnealing, iterScale, j,
														"BestSampledLogLL", bestLogLL, treeName, time));
											else
												out.println(CSV.body(m, "", "","", iterScale, j,
														"BestSampledLogLL", bestLogLL, treeName, time));
										}
										if (goldut == null)
										{
											LogInfo.logsForce("Computing gold tree using " + refMethod);
											goldut = refMethod.doIt(this, refIterScaling, goldut, treeName).getConsensus();
											computedRefTrees.mkdir();
											File current = new File(computedRefTrees, "computedRefTree_"  +f.getName().replaceAll("[.]msf","") + ".newick");
											IO.writeToDisk(current, goldut.toNewick());
										}
										for (TreeMetric tm : TreeEvaluator.coreTreeMetrics)
										{
											final double value = tm.score(inferred, goldut);
											stats.addValue(tm.toString(), value);
											if (m == InferenceMethod.ANNEALING)
												out.println(CSV.body(m, adaptiveTempDiff, adaptiveType,
														this.nAnnealing, iterScale, j, tm, value,
														treeName, time));
											else
												out.println(CSV.body(m, "", "", "", iterScale, j, tm, value,
														treeName, time));
										}

										LogInfo.end_track();
										//          ReportProgress.divisionCompleted();

										LogInfo.track("Score for current block of repeats (Method="+m + ",IterScale=" + iterScale + ",TreeName=" + treeName + ")");
										for (TreeMetric tm : TreeEvaluator.coreTreeMetrics)
											LogInfo.logsForce("Current " + tm + ":" + stats.median(tm.toString()));
										LogInfo.end_track();										
										out.flush();
										//        ReportProgress.divisionCompleted();

										if ((l < nRun - 1) && (this.nAnnealing > nIter))
										{
											nIter = this.nAnnealing;			
											nMrBayesIter=Math.max(100000, (int) (nIter*iterScalings.get(i)));
										}
										if (l == nRun - 2) {
											nAnnealing = nIter;											
											adaptiveTempDiff = false;
											adaptiveType = 0;
										}
										LogInfo.end_track();
									}									
									adaptiveTempDiff = adaptiveTempDiff0;
								}
							}
							LogInfo.end_track();

							//      ReportProgress.divisionCompleted();

		}
		//		String mrBayesFolder=Execution.getFile("temp-mrbayes");		
		//	LogInfo.logs("grep Mean: "+mrBayesFolder+"/time=*/mrbayes-stdout |awk {'print $3'}");
		//	String str =IO.call("bash -s", "grep Mean: "+mrBayesFolder+"/time=*/mrbayes-stdout |awk {'print $3'}");        
		//	logZout.append(str);
		//LogInfo.logs(str);

		out.close();
		logZout.close();
	}


	public static double numericalIntegratedMarginalLikelihood(Random rand, UnrootedTreeState ncts, double rate, int K)
	{		
		Gamma exponentialPrior = Gamma.exponential(rate);
		double[] result=new double[K];
		for(int i=0;i<K;i++)					
			result[i]=ncts.copyAndChange(AnnealingKernel.generate(rand, exponentialPrior,ncts.getUnrootedTree().leaves())).logLikelihood();
		double max = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < K; i++)
			max = Math.max(max, result[i]);
		for(int i = 0; i < K; i++)
			result[i] = Math.exp(result[i]-max);	   
		System.out.println(result);
		double finalresult=0;
		for(int i=0;i<K;i++)finalresult+=result[i];
		return  Math.log(finalresult)-Math.log(K)+max;
	}

	public static enum InferenceMethod
	{
		MCMC {
			@Override
			public TreeDistancesProcessor doIt(SMCSamplerExperiments instance, double iterScale, UnrootedTree goldut, String treeName)
			{
				PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
				samplerMain.alignmentInputFile = instance.data;
				samplerMain.st = instance.sequenceType;
				PhyloSampler._defaultPhyloSamplerOptions.nIteration = (int) (iterScale * instance.nThousandIters * 1000);
				samplerMain.run();
				return  samplerMain.tdp;
			}
		},		 
		MB {
			@Override
			public TreeDistancesProcessor doIt(SMCSamplerExperiments instance,
					double iterScale,  UnrootedTree goldut, String treeName)
			{
				MrBayes mb = MrBayes.instance;
				mb.nChains = 1;
				mb.seed = mainRand.nextInt();
				mb.nMCMCIters = (int) (iterScale * instance.nThousandIters * 1000);
				if(mb.fixGTRGammaPara)
				{
					mb.alpha=instance.generator.alpha;
					mb.subsRates=instance.generator.subsRates;
					mb.stationaryDistribution=instance.generator.stationaryDistribution;
				}
				mb.computeSamples(MSAParser.parseMSA(instance.data), instance.sequenceType);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				mb.processMrBayesTrees(tdp,1);
				mb.seed = mainRand.nextInt();
				mb.nMCMCIters = (int) (iterScale * instance.nThousandIters * 1000);
				String marginalLike= mb.computeMarginalLike(MSAParser.parseMSA(instance.data), instance.sequenceType);
				//				mb.cleanUpMrBayesOutput();
				instance.logZout.println(CSV.body(treeName,"MB", "NA",
						marginalLike));
				instance.logZout.flush();
				return tdp;
			}
		},
		ANNEALING {

			@Override
			public TreeDistancesProcessor doIt(SMCSamplerExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName)
			{
				// String resultFolder = Execution.getFile("results");
				// PrintWriter out = IOUtils.openOutEasy(new File(new File(
				// resultFolder), treeName + "logZEst.csv"));
				// out.println(CSV.header("logZ", "varLogZ"));
				SMCSampler<UnrootedTreeState> pc = new SMCSampler<UnrootedTreeState>();
				// ParticleFilter<UnrootedTreeState> pc = new
				// ParticleFilter<UnrootedTreeState>();
				pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				// pc.N = (int) iterScale;
				pc.setEss(pc.N);
				pc.nThreads = instance.nThreads;
				pc.resampleLastRound = true;
				pc.rand = instance.mainRand;
				pc.verbose = instance.verbose;
				pc.setProcessSchedule(new PhyloPFSchedule());
				// pc.resamplingStrategy =
				// ParticleFilterSMCSampler.ResamplingStrategy.ESS;
				pc.resamplingStrategy = instance.resamplingStrategy;
				pc.adaptiveTempDiff = instance.adaptiveTempDiff;
				pc.alpha = instance.alphaSMCSampler;
				pc.essRatioThreshold = instance.essRatioThreshold;
				pc.adaptiveType = instance.adaptiveType;
				String filename ="essTempDiffDeterministic.csv";
				if (instance.adaptiveTempDiff)
					filename = "essTempDiffAdaptive" + instance.adaptiveType
					+ ".csv";
				pc.smcSamplerOut = IOUtils.openOutEasy(new File(
						instance.output, filename));
				List<Taxon> leaves = MSAParser.parseMSA(instance.data).taxa();
				// RootedTree initTree = TreeGenerators.sampleCoalescent(
				// instance.mainRand, leaves, false);
				RootedTree initTree = TreeGenerators.sampleExpNonclock(
						instance.mainRand, leaves.size(), 10.0);
				// .sampleCoalescent(
				// instance.mainRand, leaves, false);
				Gamma exponentialPrior = Gamma.exponential(10.0);
				// StandardNonClockPriorDensity priorDensity = new
				// StandardNonClockPriorDensity(
				// Distrib < Double > branchDistribution);
				StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
						exponentialPrior);
				Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);		        
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
				UnrootedTreeState ncts = UnrootedTreeState.initFastState(
						UnrootedTree.fromRooted(initTree), dataset, ctmc,
						priorDensity);			
				if(dataset.observations().size()<=5)
					instance.marginalLogLike=numericalIntegratedMarginalLikelihood(instance.mainRand, ncts, 10.0, instance.nNumericalIntegration);
				System.out.println();
				System.out.println("logPrior: "+ncts.getLogPrior());
				ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
				//proposalOptions.useGlobalMultiplicativeBranchProposal=false;
				proposalOptions.useSubtreePruningRegraftingProposal=false;
				proposalOptions.multiplicativeBranchProposalScaling=1.2;
				LinkedList<ProposalDistribution> proposalDistributions = new LinkedList<ProposalDistribution>();
				// ParticleKernel<UnrootedTreeState> ppk
				AnnealingKernel ppk = new AnnealingKernel(ncts, 1.0/instance.nAnnealing, proposalDistributions, proposalOptions);				
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				pc.sample(ppk, tdp);
				String methodname=instance.adaptiveTempDiff?"Adaptive":"Deterministic";				
				instance.logZout.println(CSV.body(treeName,methodname,instance.marginalLogLike,
						pc.estimateNormalizer()));
				instance.logZout.flush();
				LogInfo.track("Estimation of log(Z) ");
				LogInfo.logsForce("Estimate of log(Z): "
						+ pc.estimateNormalizer()
						+ "; Estimate of variance of log(Z): "
						+ pc.estimateNormalizerVariance());
				LogInfo.end_track();
				instance.nAnnealing = ppk.getCurrentIter();
				return tdp;
			}
		};
		abstract  TreeDistancesProcessor doIt(SMCSamplerExperiments instance, double iterScale, UnrootedTree goldut, String treeName);
	}


	@SuppressWarnings("unchecked")
	public static <T> double dist(Counter<T> c1, Counter<T> c2)
	{
		double sum = 0;
		for (T t : union(c1.keySet(), c2.keySet()))
			sum += Math.abs(c1.getCount(t) - c2.getCount(t));
		return sum;
	}

	public static  double dist(TreeDistancesProcessor p1, TreeDistancesProcessor p2)
	{
		return dist(p1.getUnrootedCladesPosterior(), p2.getUnrootedCladesPosterior());
	}

	@Override
	public void run()
	{

		if (methods.size() != iterScalings.size())
			throw new RuntimeException("Number of methods and scaling iters should match");

		if(useDataGenerator)
		{
			LogInfo.logsForce("Generating data...");
			LogInfo.forceSilent = false;
			generator.run();
		}else{

			LogInfo.forceSilent = false;
			File dataDir = new File(Execution.getFile(dataDirName));
			LogInfo.logsForce("Copying data to "+dataDir);
			dataDir.mkdir();			
			String str =IO.call("/bin/cp " +  dataFile +" "+ dataDir+"/"+dataFile.getName());        
			LogInfo.logs(str);
		}

		//  LogInfo.forceSilent = false;
		output = new File(Execution.getFile("results"));
		output.mkdir();
		treeComparison();

	}

}
