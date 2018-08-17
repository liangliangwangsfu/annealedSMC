package evolmodel;

import static nuts.io.IO.writeToDisk;
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
import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import pty.pmcmc.PhyloPFSchedule;
import pty.smc.LazyParticleFilter;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.PriorPriorKernel;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.models.CTMC;
import smcsampler.AnnealingKernel;
import smcsampler.LinkedImportanceSampling;
import smcsampler.MrBayes;
import smcsampler.SMCSampler;
import smcsampler.SSreverse;
import smcsampler.SteppingStone;
import smcsampler.phyloMCMC2;
import smcsampler.readCSV;
import smcsampler.subSamplingKernel;
import smcsampler.subSamplingSMCsampler;
import conifer.trees.StandardNonClockPriorDensity;
import ev.ex.PhyloSamplerMain;
import ev.ex.TreeGenerators;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import ev.to.NJ;
import fig.basic.IOUtils;
import fig.basic.LogInfo;

import fig.basic.Option;
import fig.exec.Execution;
import fig.prob.Gamma;
import goblin.Taxon;

public class ModelComparisonExperiments implements Runnable
{
	@Option public boolean resampleRoot = false;
	@Option public boolean useLIS = false;
	@Option public boolean usenewSS = false;
	@Option public boolean useRevSS = false;
	@Option public boolean usenewDSMC = false;
	@Option public boolean useRef = false;
	@Option public int ntempSS = 10;
	@Option public int mcmcfac = 1;
	@Option public InferenceMethod refMethod = InferenceMethod.MB;
	@Option public static EvolutionModel evolModel = EvolutionModel.K2P;  
	@Option public double nThousandIters = 10;
	@Option public ArrayList<InferenceMethod> methods = list(Arrays.asList(InferenceMethod.ANNEALING,InferenceMethod.MB));
	@Option public ArrayList<Double> iterScalings = list(Arrays.asList(1.0));
	@Option public int refIterScaling = 100;
	@Option public int repPerDataPt = 1;    
	@Option public double pmcmcSMCExpMix = 2.0/3.0;
	@Option public double treeRatePrior=1.0;
	public static ev.ex.DataGenerator.DataGeneratorMain generator = new ev.ex.DataGenerator.DataGeneratorMain();
	public static PhyloSamplerMain samplerMain = new PhyloSamplerMain();
	public static phyloMCMC2 samplerMain2 = new phyloMCMC2();
	public static LinkedImportanceSampling LISMain = new LinkedImportanceSampling();
	public static SteppingStone ssMain = new SteppingStone();
	public static SSreverse RevssMain = new SSreverse();
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
	@Option public boolean useDataGenerator=true; 	
	@Option public File dataFile=null;
	@Option public File  refTree=null;
	@Option public String dataDirName = "output";		
	@Option public SequenceType sequenceType=SequenceType.RNA; 
	@Option public boolean isPMCMC4clock=true;
	@Option public boolean  useNJinfo=false;	
	@Option public boolean  useTopologyProcessor=false; 
	@Option public boolean  saveTreesFromPMCMC=false;
	@Option public String nameOfAllTrees="allTrees.trees";
	@Option public double csmc_trans2tranv=2.0;
	@Option public double smcmcmcMix = 0.5;
	@Option public boolean betterStartVal=true;
	@Option public SMCSampler.ResamplingStrategy resamplingStrategy=SMCSampler.ResamplingStrategy.ESS;
	@Option public subSamplingSMCsampler.ResamplingStrategy resamplingStrategy2=subSamplingSMCsampler.ResamplingStrategy.ESS;
	@Option
	public int nCSMC = 10;
	@Option
	public int nUCSMC = 10;
	@Option
	public boolean adaptiveTempDiff = false;
	@Option
	public boolean runDSMCusingadaptiveTemp = true;
	@Option
	public int adaptiveType = 0;

	@Option
	public boolean useCESS=true;

	@Option
	public double alphaSMCSampler = 0.95;

	@Option
	public double essRatioThreshold = 0.7;
	@Option public int nNumericalIntegration = 1000;
	
	@Option public int nSitesPerIndex = 10;

	@Option
	public int nAnnealing = 5000;
	@Option
	public int nSubsampling = 5000;
	public  File data = null;
	private File output = null;
	private RootedTree goldrt;
	private PrintWriter logZout = null;
	private double marginalLogLike=0;

	public static void main(String[] args)
	{
		IO.run(args, new ModelComparisonExperiments(), 
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
			CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),csmc_trans2tranv); //TODO:to be updated when the parameters are estimated 

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
//									if (m==InferenceMethod.ANNEALING && adaptiveTempDiff && runDSMCusingadaptiveTemp == true)
//										nRun = 2;
									int nIter = Integer.MIN_VALUE;
									for (int l = 0; l < nRun; l++) {
										double iterScale = iterScalings.get(i);
										if((m == InferenceMethod.MB || m == InferenceMethod.MCMC|| m == InferenceMethod.MCMC2)&& nMrBayesIter>0) iterScale=nMrBayesIter;
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
										if(m==InferenceMethod.ANNEALING && l==0) {
											//IO.writeToDisk(new File(output, "consensus_annealing_"+treeName.replaceAll("[.]msf$",".newick")), inferred.toNewick());
											IO.writeToDisk(new File(output, "consensus_annealing.newick"), inferred.toNewick());
										}
										IO.writeToDisk(new File(output, "consensus_"+treeName.replaceAll("[.]msf$",".newick")), inferred.toNewick());
										{
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
											nMrBayesIter=Math.max(1000, (int) (nIter*iterScalings.get(i)));
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
		CSMC {

			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut,
					String treeName) {
				ParticleFilterOptions options = new ParticleFilterOptions();

				// ParticleFilter<PartialCoalescentState> pc = new
				// ParticleFilter<PartialCoalescentState>();
				options.nParticles = (int) (iterScale * instance.nThousandIters * 1000);
				options.nThreads = instance.nThreads;
				options.resampleLastRound = true;
				options.parallelizeFinalParticleProcessing = true;
				options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
				options.maxNUniqueParticles = instance.maxNUniqueParticles;
				options.rand = instance.mainRand;
				options.verbose = instance.verbose;
				// options.maxNGrow = 0;

				Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);				
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), instance.csmc_trans2tranv);
				PartialCoalescentState init = PartialCoalescentState.initFastState(dataset, ctmc);
				LazyParticleKernel pk2 = new PriorPriorKernel(init);
				LazyParticleFilter<PartialCoalescentState> lpf = new LazyParticleFilter<PartialCoalescentState>(pk2,
						options);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();

				if (instance.useTopologyProcessor) {
					TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
					final double zHat = lpf.sample(tdp, trTopo);

					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: " + urtCounter.keySet().size());
					for (UnrootedTree urt : urtCounter.keySet()) {
						LogInfo.logsForce(urt);
						LogInfo.logsForce(urtCounter.getCount(urt));
					}
				} else {
					final double zHat = lpf.sample(tdp);
					// LogInfo.logsForce("Norm:" + zHat);
				}
				return tdp;
			}
		},
		CSMCNonClock {

			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut,
					String treeName) {
				ParticleFilterOptions options = new ParticleFilterOptions();
				// ParticleFilter<PartialCoalescentState> pc = new
				// ParticleFilter<PartialCoalescentState>();
				options.nParticles = (int) (iterScale * instance.nThousandIters * 1000);
				options.nThreads = instance.nThreads;
				options.resampleLastRound = true;
				options.parallelizeFinalParticleProcessing = true;
				options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
				options.maxNUniqueParticles = instance.maxNUniqueParticles;
				options.rand = instance.mainRand;
				options.verbose = instance.verbose;
				// options.maxNGrow = 0;
				Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(), instance.csmc_trans2tranv);
				PartialCoalescentState init = PartialCoalescentState.initFastState(instance.resampleRoot, dataset, ctmc,
						false);
				LazyParticleKernel pk2 = new NCPriorPriorKernel(init);
				LazyParticleFilter<PartialCoalescentState> lpf = new LazyParticleFilter<PartialCoalescentState>(pk2,
						options);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				final double zHat = lpf.sample(tdp);
				
				String methodname = "CSMCNonClock";				
				instance.logZout.println(CSV.body(treeName,methodname,zHat));
				instance.logZout.flush();
				
				LogInfo.logsForce("Norm:" + zHat);
				return tdp;
			}
		},
		MCMC {
			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut, String treeName)
			{

				PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
				samplerMain.alignmentInputFile = instance.data;
				samplerMain.st = instance.sequenceType;
				int Ntemperature = 4;
				double alpha = 1.0/3.0;
				PhyloSampler._defaultPhyloSamplerOptions.nIteration = (int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.mcmcfac));
				//PhyloSampler._defaultPhyloSamplerOptions.nIteration = 10000;
				//double[] temperatureSchedule = new double[]{1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0};	
				samplerMain.setTemperatureSchedule(Ntemperature, alpha);
				samplerMain.computeLogZUsingSteppingStone = false;
				samplerMain.setnSamplesEachChain(PhyloSampler._defaultPhyloSamplerOptions.nIteration);
				samplerMain.setLogZ(0.0);
				samplerMain.run();			
				
				instance.logZout.println(CSV.body(treeName,"MCMC", "NA",
						samplerMain.getLogZ()));
				instance.logZout.flush();
				
				
				if(instance.useLIS == true) {
					LinkedImportanceSampling._defaultPhyloSamplerOptions.rand = mainRand;
					LISMain.alignmentInputFile = instance.data;
					LISMain.st = instance.sequenceType;
					LISMain.nSamplesEachChain = ((int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.ntempSS)));
					//LISMain.nSamplesEachChain = 10000;
					LISMain.setnSamplesEachChain(LISMain.nSamplesEachChain);
					LISMain.nChains = instance.ntempSS;
					LISMain.alpha = 1.0/3.0;
					LISMain.csmc_trans2tranv = instance.csmc_trans2tranv;
					LISMain.run();	
					instance.logZout.println(CSV.body(treeName,"LIS", "NA",
							LISMain.getNormalizer()));
					instance.logZout.flush();	
				}
				
				if(instance.useRevSS == true)  {
					SSreverse._defaultPhyloSamplerOptions.rand = mainRand;
					RevssMain.alignmentInputFile = instance.data;
					RevssMain.st = instance.sequenceType;
					RevssMain.nSamplesEachChain = ((int) (iterScale * instance.nThousandIters * 1000/(1.0*(instance.ntempSS+1))));
					RevssMain.setnSamplesEachChain(RevssMain.nSamplesEachChain);
					RevssMain.nChains = instance.ntempSS;
					RevssMain.alpha = 1.0/3.0;
					RevssMain.csmc_trans2tranv = instance.csmc_trans2tranv;
					RevssMain.run();	
					instance.logZout.println(CSV.body(treeName,"RevSS", "NA",
							RevssMain.getNormalizer()));
					instance.logZout.flush();	
				}
				
				if(instance.usenewSS == true){
					SteppingStone._defaultPhyloSamplerOptions.rand = mainRand;
					ssMain.alignmentInputFile = instance.data;
					ssMain.st = instance.sequenceType;
					ssMain.nSamplesEachChain = ((int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.ntempSS)));
					ssMain.setnSamplesEachChain(ssMain.nSamplesEachChain);
					ssMain.nChains = instance.ntempSS;
					ssMain.alpha = 1.0/3.0;
					ssMain.csmc_trans2tranv = instance.csmc_trans2tranv;
					ssMain.run();	
					instance.logZout.println(CSV.body(treeName,"SS", "NA",
							ssMain.getNormalizer()));
					instance.logZout.flush();	
				}
									
				return  samplerMain.tdp;
			}
		},	
		MCMC2 {
			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut, String treeName)
			{

				//PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
				PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
				samplerMain2.alignmentInputFile = instance.data;
				samplerMain2.st = instance.sequenceType;
				//String newname = instance.output+"/"+"consensus_annealing_" +treeName.replaceAll("[.]msf$",".newick");
				String newname = instance.output+"/"+ "consensus_annealing.newick";
				samplerMain2.initTree =  UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(new File(newname))));
				PhyloSampler._defaultPhyloSamplerOptions.nIteration = (int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.mcmcfac));
				samplerMain2.run();			
				
				return  samplerMain2.tdp;
			}
		},
		MB {
			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance,
					double iterScale,  UnrootedTree goldut, String treeName)
			{
				MrBayes mb = MrBayes.instance;
				mb.nChains = 1;
				mb.seed = mainRand.nextInt();
				mb.nMCMCIters = (int) (iterScale * instance.nThousandIters * 1000);
				List<Taxon> leaves = MSAParser.parseMSA(instance.data).taxa();
				UnrootedTree initTree = initTree(new Random(5),  leaves);
				String outName="startTree.newick";
				writeToDisk(new File(instance.output, outName), initTree.toNewick());
				String cmdStr="cat " +outName + "  | sed 's/internal_[0-9]*//g' > " + "start-tree.newick";                    
//				LogInfo.logs(cmdStr);
                IO.call("bash -s",cmdStr,instance.output);          
//				mb.setStartTree(IO.f2s(new File(instance.output,"start-tree.newick")));
				
				if(leaves.size()<4) mb.useNNI=false;
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
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance,
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
				pc.setUseCESS(instance.useCESS);
				if(pc.adaptiveTempDiff == false && instance.usenewDSMC == false) {
					 String str =instance.output+"/"+"essTempDiffAdaptive00.csv";
					 pc.setDeterministicTemperatureDifference(readCSV.DeterministicTem(str));
				}
				if(instance.usenewDSMC == true) {
					int T = instance.nAnnealing;
					List<Double> DeterministicTD = new ArrayList<Double>();
					DeterministicTD.add(0.0);
					for(int t = 1; t < T+1; t++) {
						DeterministicTD.add(Math.pow((t*1.0)/(T*1.0), 3)-Math.pow(((t-1)*1.0)/(T*1.0), 3));
					}
					pc.setDeterministicTemperatureDifference(DeterministicTD);
				}
				String filename ="essTempDiffDeterministic.csv";
				String filename2 ="essTempDiffDeterministic00.csv";
				if (instance.adaptiveTempDiff) {
					filename = "essTempDiffAdaptive" + instance.adaptiveType
							+ ".csv";
					filename2 ="essTempDiffAdaptive00.csv";
				}
					
				pc.smcSamplerOut = IOUtils.openOutEasy(new File(
						instance.output, filename));
				pc.smcSamplerOut2 = IOUtils.openOutEasy(new File(
						instance.output, filename2));
				List<Taxon> leaves = MSAParser.parseMSA(instance.data).taxa();
				UnrootedTree initTree = initTree(new Random(3), leaves);
				Gamma exponentialPrior = Gamma.exponential(10.0);
				StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
						exponentialPrior);
				Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);		        
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),instance.csmc_trans2tranv);
				UnrootedTreeState ncts = UnrootedTreeState.initFastState(initTree, dataset, ctmc, priorDensity);			
				if(dataset.observations().size()<=5)
					instance.marginalLogLike=numericalIntegratedMarginalLikelihood(instance.mainRand, ncts, 10.0, instance.nNumericalIntegration);

				LogInfo.track("log prior and loglikelihood of the initial tree: ");
				LogInfo.logsForce(ncts.getLogPrior()+" "+ncts.getLogLikelihood());
				LogInfo.end_track();

				ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
				if(MSAParser.parseMSA(instance.data).nTaxa()<4)  proposalOptions.useStochasticNearestNeighborInterchangeProposal=false;
				else
					proposalOptions.useStochasticNearestNeighborInterchangeProposal=true;				
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
		},
		ANNEALINGEvolPara {

			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName)
			{
				// String resultFolder = Execution.getFile("results");
				// PrintWriter out = IOUtils.openOutEasy(new File(new File(
				// resultFolder), treeName + "logZEst.csv"));
				// out.println(CSV.header("logZ", "varLogZ"));
				SMCSampler<UnrootedTreeEvolParameterState> pc = new SMCSampler<UnrootedTreeEvolParameterState>();
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
				pc.setUseCESS(instance.useCESS);
				if(pc.adaptiveTempDiff == false && instance.usenewDSMC == false) {
					 String str =instance.output+"/"+"essTempDiffAdaptive00.csv";
					 pc.setDeterministicTemperatureDifference(readCSV.DeterministicTem(str));
				}
				if(instance.usenewDSMC == true) {
					int T = instance.nAnnealing;
					List<Double> DeterministicTD = new ArrayList<Double>();
					DeterministicTD.add(0.0);
					for(int t = 1; t < T+1; t++) {
						DeterministicTD.add(Math.pow((t*1.0)/(T*1.0), 3)-Math.pow(((t-1)*1.0)/(T*1.0), 3));
					}
					pc.setDeterministicTemperatureDifference(DeterministicTD);
				}
				String filename ="essTempDiffDeterministic.csv";
				String filename2 ="essTempDiffDeterministic00.csv";
				if (instance.adaptiveTempDiff) {
					filename = "essTempDiffAdaptive" + instance.adaptiveType
							+ ".csv";
					filename2 ="essTempDiffAdaptive00.csv";
				}
					
				pc.smcSamplerOut = IOUtils.openOutEasy(new File(
						instance.output, filename));
				pc.smcSamplerOut2 = IOUtils.openOutEasy(new File(
						instance.output, filename2));
				List<Taxon> leaves = MSAParser.parseMSA(instance.data).taxa();
				UnrootedTree initTree = initTree(new Random(3), leaves);
				Gamma exponentialPrior = Gamma.exponential(10.0);
				StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
						exponentialPrior);
				Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);		        				
			//	EvolutionParameters evolPara = new EvolutionParameters.K2P(instance.csmc_trans2tranv);								
				EvolutionParameters evolPara = new EvolutionParameters.GTR(new double[]{0.26, 0.18, 0.17, 0.15, 0.11, 0.13, 0.25, 0.25, 0.25, 0.25});
				
				CTMC ctmc = evolModel.instantiateCTMC(evolPara, dataset.nSites());				
				UnrootedTreeState ncts = UnrootedTreeState.initFastState(initTree, dataset, ctmc, priorDensity);			
				if(dataset.observations().size()<=5)
					instance.marginalLogLike=numericalIntegratedMarginalLikelihood(instance.mainRand, ncts, 10.0, instance.nNumericalIntegration);

				LogInfo.track("log prior and loglikelihood of the initial tree: ");
				LogInfo.logsForce(ncts.getLogPrior()+" "+ncts.getLogLikelihood());
				LogInfo.end_track();

				ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
				if(MSAParser.parseMSA(instance.data).nTaxa()<4)  proposalOptions.useStochasticNearestNeighborInterchangeProposal=false;
				else
					proposalOptions.useStochasticNearestNeighborInterchangeProposal=true;		
				LinkedList<ProposalDistribution> proposalDistributions = new LinkedList<ProposalDistribution>();
				
				EvolutionParameterProposalDistribution.Options evolProposalOptions = EvolutionParameterProposalDistribution.Util._defaultProposalDistributionOptions;
				evolProposalOptions.useGTRProposal = true;
				evolProposalOptions.useK2PProposal = false;
				LinkedList<EvolutionParameterProposalDistribution> evolProposalDistributions = new LinkedList<EvolutionParameterProposalDistribution>();
								
				AnnealingKernelTreeEvolPara ppk = new AnnealingKernelTreeEvolPara(dataset, priorDensity, new UnrootedTreeEvolParameterState(ncts, evolPara), 1.0/instance.nAnnealing, proposalDistributions, proposalOptions, evolProposalDistributions, evolProposalOptions);
				ppk.setEvolModel(evolModel);
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
		},
		SUBSAMPLING {
			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName)
			{
				// String resultFolder = Execution.getFile("results");
				// PrintWriter out = IOUtils.openOutEasy(new File(new File(
				// resultFolder), treeName + "logZEst.csv"));
				// out.println(CSV.header("logZ", "varLogZ"));
				subSamplingSMCsampler<UnrootedTreeState> pc = new subSamplingSMCsampler<UnrootedTreeState>();
				// ParticleFilter<UnrootedTreeState> pc = new
				// ParticleFilter<UnrootedTreeState>();
				pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				// pc.N = (int) iterScale;
				pc.setEss(pc.N);
				pc.nThreads = instance.nThreads;
				pc.resampleLastRound = false;
				pc.rand = instance.mainRand;
				pc.verbose = instance.verbose;
				pc.setProcessSchedule(new PhyloPFSchedule());
				// pc.resamplingStrategy =
				// ParticleFilterSMCSampler.ResamplingStrategy.ESS;
				pc.resamplingStrategy = instance.resamplingStrategy2;
				//pc.adaptiveTempDiff = instance.adaptiveTempDiff;
				//pc.alpha = instance.alphaSMCSampler;
				pc.essRatioThreshold = instance.essRatioThreshold;
				//pc.adaptiveType = instance.adaptiveType;
				pc.setUseCESS(instance.useCESS);
//				if(pc.adaptiveTempDiff == false && instance.usenewDSMC == false) {
//					 String str =instance.output+"/"+"essTempDiffAdaptive00.csv";
//					 pc.setDeterministicTemperatureDifference(readCSV.DeterministicTem(str));
//				}
//				if(instance.usenewDSMC == true) {

				
				//}
//				String filename ="essTempDiffDeterministic.csv";
//				String filename2 ="essTempDiffDeterministic00.csv";
//				if (instance.adaptiveTempDiff) {
//					filename = "essTempDiffAdaptive" + instance.adaptiveType
//							+ ".csv";
//					filename2 ="essTempDiffAdaptive00.csv";
//				}
					
//				pc.smcSamplerOut = IOUtils.openOutEasy(new File(
//						instance.output, filename));
//				pc.smcSamplerOut2 = IOUtils.openOutEasy(new File(
//						instance.output, filename2));
				int T = instance.nSubsampling;
				List<Double> originalTemp = new ArrayList<Double>();
				if(instance.adaptiveTempDiff == true) {
				   //String str =instance.output+"/"+"essTempDiffAdaptive00.csv";
				   String str ="/Users/oudomame/Desktop/essTempDiffAdaptive3.csv";
				   List<Double> temperatureDifferenceFromAnnealing = readCSV.DeterministicTem(str);
				   T = temperatureDifferenceFromAnnealing.size()-1;			   
				   originalTemp.add(temperatureDifferenceFromAnnealing.get(1));
					for(int t = 1; t < T; t++) {
						originalTemp.add(originalTemp.get(t-1) + temperatureDifferenceFromAnnealing.get(t+1));
					}		   
			    }else {
					for(int t = 0; t < T; t++) {
						originalTemp.add(Math.pow((1.0*t+1)/(T*1.0), 3));
					}		    	
			    }
				
				
				List<Taxon> leaves = MSAParser.parseMSA(instance.data).taxa();
				UnrootedTree initTree = initTree(new Random(3), leaves);
				Gamma exponentialPrior = Gamma.exponential(10.0);
				StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
						exponentialPrior);
				Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);	
				int nSitesPerIndex = instance.nSitesPerIndex;
				int nFakeSites = 0;
				int nSites = dataset.nSites();
				if((nSites/nSitesPerIndex)*nSitesPerIndex == nSites) {
					nFakeSites = nSites/nSitesPerIndex;
				}else {
					nFakeSites = nSites/nSitesPerIndex + 1;
				}
				
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),instance.csmc_trans2tranv);
				UnrootedTreeState ncts = UnrootedTreeState.initFastState(initTree, dataset, ctmc, priorDensity);			
//				if(dataset.observations().size()<=5)
//					instance.marginalLogLike=numericalIntegratedMarginalLikelihood(instance.mainRand, ncts, 10.0, instance.nNumericalIntegration);

				LogInfo.track("log prior and loglikelihood of the initial tree: ");
				LogInfo.logsForce(ncts.getLogPrior()+" "+ncts.getLogLikelihood());
				LogInfo.end_track();

				ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
				if(MSAParser.parseMSA(instance.data).nTaxa()<4)  proposalOptions.useStochasticNearestNeighborInterchangeProposal=false;
				else
					proposalOptions.useStochasticNearestNeighborInterchangeProposal=true;				
				LinkedList<ProposalDistribution> proposalDistributions = new LinkedList<ProposalDistribution>();
				// ParticleKernel<UnrootedTreeState> ppk
				subSamplingKernel ppk = new subSamplingKernel(ncts, T, proposalDistributions, proposalOptions);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				pc.setData(dataset);
				pc.setnSites(dataset.nSites());
				pc.setnFakeSites(nFakeSites);
				pc.setnSitesPerIndex(nSitesPerIndex);
				//pc.setNewTemperature2(originalTemp, dataset.nSites());
				pc.setNewTemperature2(originalTemp, nFakeSites);
				//pc.setNewTemperatureDiff(dataset.nSites());
				pc.sample(ppk, tdp);
				String methodname = "Subsampling";				
				instance.logZout.println(CSV.body(treeName,methodname,instance.marginalLogLike,
						pc.estimateNormalizer()));
				instance.logZout.flush();
				LogInfo.track("Estimation of log(Z) ");
				LogInfo.logsForce("Estimate of log(Z): "
						+ pc.estimateNormalizer()
						+ "; Estimate of variance of log(Z): "
						+ pc.estimateNormalizerVariance());
				LogInfo.end_track();
				//instance.nAnnealing = ppk.getCurrentIter();
				return tdp;
			}		
		};
		abstract  TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut, String treeName);
	}





	public static  UnrootedTree initTree(Random rand, List<Taxon> leaves){	
		return UnrootedTree.fromRooted(TreeGenerators.sampleExpNonclock(rand,leaves, 10.0));				
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
