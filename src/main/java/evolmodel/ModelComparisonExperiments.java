package evolmodel;

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
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.PriorPriorKernel;
import pty.smc.models.CTMC;
import smcsampler.AnnealingKernel;
import smcsampler.LinkedImportanceSampling;
import smcsampler.MrBayes;
import smcsampler.SMCSampler;
import smcsampler.SSreverse;
import smcsampler.SteppingStone;
import smcsampler.phyloMCMC2;
import smcsampler.readCSV;
import smcsampler.subSamplingSMCsampler;
import conifer.trees.StandardNonClockPriorDensity;
import dr.math.distributions.GammaDistribution;
import ev.ex.PhyloSamplerMain;
import ev.ex.TreeGenerators;
import ev.poi.processors.TreeDistancesProcessor;
import ev.to.NJ;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import fig.prob.Dirichlet;
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
	@Option public InferenceMethod refMethod = InferenceMethod.SSJC;  
	@Option public double nThousandIters = 10;
	@Option public ArrayList<InferenceMethod> methods = list(Arrays.asList(InferenceMethod.ANNEALINGJC,InferenceMethod.SSJC));
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
		logZout.println(CSV.header("treeName", "Method",  "logZ"));
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

		PrintWriter gtrGammaParameters = IOUtils.openOutEasy(new File(Execution.getFile("results"),
				"gtrGammaParameters.csv"));
		gtrGammaParameters.println(CSV.header("treeName", "pi_1", "pi_2", "pi_3", "pi_4", "r_1", "r_2", "r_3", "r_4", "r_5", "r_6", "alpha"));			
		int nMCMCIter=0;		
		for (File f : files)
		{
			this.data = f;
			final String treeName = f.getName();
			LogInfo.track("Current tree:" + treeName);						
			// evaluate the likelihood of the inferred tree
			Dataset dataset = DatasetUtils.fromAlignment(this.data, sequenceType);
			UnrootedTree goldut = 
					(generator.useGutellData ||!useDataGenerator)?
							(refTree==null?null:UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(refTree)))):
								UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(new File(f.getAbsolutePath().replaceAll("[.]msf$",".newick")))));
							File computedRefTrees = new File(Execution.getFile("computed-ref-trees"));
							for (int j = 0; j < repPerDataPt; j++)
							{
								CTMC ctmc = null;
								String treeNameCurrentRep= treeName;
								treeNameCurrentRep = treeNameCurrentRep + ".Rep"+ j;
								LogInfo.track("Repeat " + (j+1) + "/" + repPerDataPt);
								LogInfo.end_track();
								gtrGammaParameters.println(CSV.body(treeNameCurrentRep,generator.stationaryDistribution[0],generator.stationaryDistribution[1], generator.stationaryDistribution[3], generator.stationaryDistribution[3],
										generator.subsRates[0],generator.subsRates[1],generator.subsRates[2],generator.subsRates[3],generator.subsRates[4],generator.subsRates[5],generator.alpha));
								gtrGammaParameters.flush();

								for (int i = 0; i < methods.size(); i++) {
									InferenceMethod m = methods.get(i);
									EvolutionParameters evolPara = null;
									EvolutionModel evolModel = EvolutionModel.JC;																			
									if(m == InferenceMethod.SSEvolK2P || m == InferenceMethod.ANNEALINGEvolK2P) {
										evolModel = EvolutionModel.K2P;
										evolPara = new EvolutionParameters.K2P(csmc_trans2tranv); 	
									}
									if(m == InferenceMethod.SSEvolGTR || m == InferenceMethod.ANNEALINGEvolGTR)
									{
										evolModel = EvolutionModel.GTR;
										evolPara = new EvolutionParameters.GTR(new double[]{0.26, 0.18, 0.17, 0.15, 0.11, 0.13, 0.25, 0.25, 0.25, 0.25});  // GTR model																								
									}
									ctmc = evolModel.instantiateCTMC(evolPara, dataset.nSites());

									double iterScale = iterScalings.get(i);
									if((m == InferenceMethod.SSJC || m == InferenceMethod.SSEvolK2P || m == InferenceMethod.SSEvolGTR)&& nMCMCIter>0) iterScale=nMCMCIter;
									LogInfo.track("Current method:" + m + " with iterScale=" + iterScale + " (i.e. " + (iterScale * nThousandIters * 1000.0) + " iterations)");
									DescriptiveStatisticsMap<String> stats = new DescriptiveStatisticsMap<String>();

									long  time=System.currentTimeMillis();
									TreeDistancesProcessor processor = m.doIt(this, iterScale, goldut, treeNameCurrentRep);									
									time=System.currentTimeMillis()-time;
									LogInfo.forceSilent = false;
									UnrootedTree inferred = processor.getConsensus();		
									System.out.println(inferred);
									if(m==InferenceMethod.ANNEALINGJC) {
										IO.writeToDisk(new File(output, "consensus_annealing.newick"), inferred.toNewick());
									}
									IO.writeToDisk(new File(output, "consensus_"+treeName.replaceAll("[.]msf$",".newick")), inferred.toNewick());
									{
										UnrootedTreeState ncs = UnrootedTreeState.initFastState(inferred, dataset, ctmc);
										if (m == InferenceMethod.ANNEALINGJC || m == InferenceMethod.ANNEALINGEvolK2P  || m == InferenceMethod.ANNEALINGEvolGTR)
											out.println(CSV.body(m, adaptiveTempDiff, adaptiveType,
													this.nAnnealing, iterScale, j,
													"ConsensusLogLL", ncs.logLikelihood(),
													treeName, time));
										else
											out.println(CSV.body(m, "", "","", iterScale, j,
													"", "",
													treeName, time));											
									}
									{
										// best log likelihood, when available
										double bestLogLL = processor.getBestLogLikelihood();
										if (m == InferenceMethod.ANNEALINGJC  || m == InferenceMethod.ANNEALINGEvolK2P  || m == InferenceMethod.ANNEALINGEvolGTR)
											out.println(CSV.body(m, adaptiveTempDiff, adaptiveType,
													this.nAnnealing, iterScale, j,
													"BestSampledLogLL", bestLogLL, treeName, time));
										//										else
										//											out.println(CSV.body(m, "", "","", iterScale, j,
										//													"BestSampledLogLL", bestLogLL, treeName, time));
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
										if (m == InferenceMethod.ANNEALINGJC  ||  m == InferenceMethod.ANNEALINGEvolK2P  || m == InferenceMethod.ANNEALINGEvolGTR)
											out.println(CSV.body(m, adaptiveTempDiff, adaptiveType,
													this.nAnnealing, iterScale, j, tm, value,
													treeName, time));
										//else
										//	out.println(CSV.body(m, "", "", "", iterScale, j, tm, value,
										//			treeName, time));
									}

								//	LogInfo.end_track();
									//          ReportProgress.divisionCompleted();

									LogInfo.track("Score for current block of repeats (Method="+m + ",IterScale=" + iterScale + ",TreeName=" + treeName + ")");
									for (TreeMetric tm : TreeEvaluator.coreTreeMetrics)
										LogInfo.logsForce("Current " + tm + ":" + stats.median(tm.toString()));
									LogInfo.end_track();										
									out.flush();
									//        ReportProgress.divisionCompleted();																				
									nMCMCIter=Math.max(1000, (int) (this.nAnnealing*iterScalings.get(i)));
									LogInfo.end_track();
								}																	
							}
							LogInfo.end_track();
		}
		out.close();
		logZout.close();
		gtrGammaParameters.close();
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
		SSJC {
			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut, String treeName)
			{

				PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
				samplerMain.alignmentInputFile = instance.data;
				samplerMain.st = instance.sequenceType;				
				int Ntemperature = 4;
				double alpha = 1.0/3.0;
				PhyloSampler._defaultPhyloSamplerOptions.nIteration = 10; //(int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.mcmcfac));
				samplerMain.setTemperatureSchedule(Ntemperature, alpha);
				samplerMain.computeLogZUsingSteppingStone = false;
				samplerMain.setnSamplesEachChain(PhyloSampler._defaultPhyloSamplerOptions.nIteration);
				samplerMain.setLogZ(0.0);
				samplerMain.run();			
				instance.logZout.flush();

				if(instance.usenewSS == true){
					SteppingStone._defaultPhyloSamplerOptions.rand = mainRand;
					ssMain.alignmentInputFile = instance.data;
					ssMain.st = instance.sequenceType;
					ssMain.nSamplesEachChain = ((int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.ntempSS)));
					ssMain.setnSamplesEachChain(ssMain.nSamplesEachChain);
					ssMain.nChains = instance.ntempSS;
					ssMain.alpha = 1.0/3.0;
					ssMain.csmc_trans2tranv = 1.0;//instance.csmc_trans2tranv;
					ssMain.run();	
					instance.logZout.println(CSV.body(treeName,"SSJC", ssMain.getNormalizer()));
					instance.logZout.flush();	
				}

				return  samplerMain.tdp;
			}
		},	
		SSEvolK2P {
			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut, String treeName)
			{
				PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
				samplerMain.alignmentInputFile = instance.data;
				samplerMain.st = instance.sequenceType;
				int Ntemperature = 4;
				double alpha = 1.0/3.0;
				//TODO: add the evolutionary parameter estimate in the MCMC.  
				PhyloSampler._defaultPhyloSamplerOptions.nIteration = 10; // (int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.mcmcfac));
				samplerMain.setTemperatureSchedule(Ntemperature, alpha);
				samplerMain.computeLogZUsingSteppingStone = false;
				samplerMain.setnSamplesEachChain(PhyloSampler._defaultPhyloSamplerOptions.nIteration);
				samplerMain.setLogZ(0.0);
				samplerMain.run();			
				instance.logZout.flush();

				if(instance.usenewSS == true){
					SteppingStoneEvolPara ssMain = new SteppingStoneEvolPara();		
					ssMain.evolModel =  EvolutionModel.K2P;  // instance.evolModel;
					Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);
					ssMain.dataset = dataset;
					ssMain.rand = mainRand;
					SteppingStone._defaultPhyloSamplerOptions.rand = mainRand;
					ssMain.alignmentInputFile = instance.data;
					ssMain.st = instance.sequenceType;
					ssMain.nSamplesEachChain = ((int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.ntempSS)));
					ssMain.setnSamplesEachChain(ssMain.nSamplesEachChain);
					ssMain.nChains = instance.ntempSS;
					ssMain.alpha = 1.0/3.0;
					ssMain.csmc_trans2tranv = instance.csmc_trans2tranv;
					ssMain.run();	
					instance.logZout.println(CSV.body(treeName,"SSEvolK2P", ssMain.getNormalizer()));
					instance.logZout.flush();	
				}

				return  samplerMain.tdp;
			}
		},
		SSEvolGTR {
			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance, double iterScale, UnrootedTree goldut, String treeName)
			{
				PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
				samplerMain.alignmentInputFile = instance.data;
				samplerMain.st = instance.sequenceType;
				int Ntemperature = 4;
				double alpha = 1.0/3.0;
				//TODO: add the evolutionary parameter estimate in the MCMC.  
				PhyloSampler._defaultPhyloSamplerOptions.nIteration = 10; // (int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.mcmcfac));
				samplerMain.setTemperatureSchedule(Ntemperature, alpha);
				samplerMain.computeLogZUsingSteppingStone = false;
				samplerMain.setnSamplesEachChain(PhyloSampler._defaultPhyloSamplerOptions.nIteration);
				samplerMain.setLogZ(0.0);
				samplerMain.run();			
				instance.logZout.flush();
				if(instance.usenewSS == true){
					SteppingStoneEvolPara ssMain = new SteppingStoneEvolPara();		
					ssMain.evolModel =  EvolutionModel.GTR;  // instance.evolModel;
					Dataset dataset = DatasetUtils.fromAlignment(instance.data, instance.sequenceType);
					ssMain.dataset = dataset;
					ssMain.rand = mainRand;
					SteppingStone._defaultPhyloSamplerOptions.rand = mainRand;
					ssMain.alignmentInputFile = instance.data;
					ssMain.st = instance.sequenceType;
					ssMain.nSamplesEachChain = ((int) (iterScale * instance.nThousandIters * 1000/(1.0*instance.ntempSS)));
					ssMain.setnSamplesEachChain(ssMain.nSamplesEachChain);
					ssMain.nChains = instance.ntempSS;
					ssMain.alpha = 1.0/3.0;
					ssMain.csmc_trans2tranv = instance.csmc_trans2tranv;
					ssMain.run();	
					instance.logZout.println(CSV.body(treeName,"SSEvolGTR", ssMain.getNormalizer()));
					instance.logZout.flush();	
				}

				return  samplerMain.tdp;
			}
		},	
		ANNEALINGJC {

			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName)
			{
				SMCSampler<UnrootedTreeState> pc = new SMCSampler<UnrootedTreeState>();
				pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				pc.setEss(pc.N);
				pc.nThreads = instance.nThreads;
				pc.resampleLastRound = true;
				pc.rand = instance.mainRand;
				pc.verbose = instance.verbose;
				pc.setProcessSchedule(new PhyloPFSchedule());
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
				//CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),instance.csmc_trans2tranv);
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),1.0);
				UnrootedTreeState ncts = UnrootedTreeState.initFastState(initTree, dataset, ctmc, priorDensity);			

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
				//String methodname=instance.adaptiveTempDiff?"Adaptive":"Deterministic";				
				instance.logZout.println(CSV.body(treeName,"ANNEALINGJC", pc.estimateNormalizer()));
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
		ANNEALINGEvolK2P {

			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName)
			{
				SMCSampler<UnrootedTreeEvolParameterState> pc = new SMCSampler<UnrootedTreeEvolParameterState>();
				pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				pc.setEss(pc.N);
				pc.nThreads = instance.nThreads;
				pc.resampleLastRound = true;
				pc.rand = instance.mainRand;
				pc.verbose = instance.verbose;
				pc.setProcessSchedule(new PhyloPFSchedule());
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
				EvolutionParameters evolPara = new EvolutionParameters.K2P(instance.csmc_trans2tranv);
				CTMC ctmc = EvolutionModel.K2P.instantiateCTMC(evolPara, dataset.nSites());
				UnrootedTreeState ncts = UnrootedTreeState.initFastState(initTree, dataset, ctmc, priorDensity);						
				LogInfo.track("log prior and loglikelihood of the initial tree: ");
				LogInfo.logsForce(ncts.getLogPrior()+" "+ncts.getLogLikelihood());
				LogInfo.end_track();

				ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
				if(MSAParser.parseMSA(instance.data).nTaxa()<4)  proposalOptions.useStochasticNearestNeighborInterchangeProposal=false;
				else
					proposalOptions.useStochasticNearestNeighborInterchangeProposal=true;		
				LinkedList<ProposalDistribution> proposalDistributions = new LinkedList<ProposalDistribution>();

				EvolutionParameterProposalDistribution.Options evolProposalOptions = EvolutionParameterProposalDistribution.Util._defaultProposalDistributionOptions;
				evolProposalOptions.useGTRProposal = false;
				evolProposalOptions.useK2PProposal = true;
				LinkedList<EvolutionParameterProposalDistribution> evolProposalDistributions = new LinkedList<EvolutionParameterProposalDistribution>();								
				AnnealingKernelTreeEvolPara ppk = new AnnealingKernelTreeEvolPara(dataset, priorDensity, new UnrootedTreeEvolParameterState(ncts, evolPara), 1.0/instance.nAnnealing, proposalDistributions, proposalOptions, evolProposalDistributions, evolProposalOptions);
				ppk.setEvolModel(EvolutionModel.K2P);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				pc.sample(ppk, tdp);
				//				String methodname=instance.adaptiveTempDiff?"Adaptive":"Deterministic";				
				instance.logZout.println(CSV.body(treeName,"ANNEALINGEvolK2P", pc.estimateNormalizer()));
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
		ANNEALINGEvolGTR {

			@Override
			public TreeDistancesProcessor doIt(ModelComparisonExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName)
			{
				SMCSampler<UnrootedTreeEvolParameterState> pc = new SMCSampler<UnrootedTreeEvolParameterState>();
				pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				pc.setEss(pc.N);
				pc.nThreads = instance.nThreads;
				pc.resampleLastRound = true;
				pc.rand = instance.mainRand;
				pc.verbose = instance.verbose;
				pc.setProcessSchedule(new PhyloPFSchedule());
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
				EvolutionParameters evolPara = new EvolutionParameters.GTR(new double[]{0.26, 0.18, 0.17, 0.15, 0.11, 0.13, 0.25, 0.25, 0.25, 0.25});  // GTR model				
				CTMC ctmc = EvolutionModel.GTR.instantiateCTMC(evolPara, dataset.nSites());
				UnrootedTreeState ncts = UnrootedTreeState.initFastState(initTree, dataset, ctmc, priorDensity);			
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
				ppk.setEvolModel(EvolutionModel.GTR);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				pc.sample(ppk, tdp);
				//				String methodname=instance.adaptiveTempDiff?"Adaptive":"Deterministic";				
				instance.logZout.println(CSV.body(treeName,"ANNEALINGEvolGTR", pc.estimateNormalizer()));
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

		output = new File(Execution.getFile("results"));
		output.mkdir();

		if(useDataGenerator)
		{
			LogInfo.logsForce("Generating data...");
			LogInfo.forceSilent = false;
			generator.stationaryDistribution = Dirichlet.sample(generator.rand, new double[] {100,100,100,100}); 
			generator.subsRates = Dirichlet.sample(generator.rand, new double[] {100,100,100,100,100,100}); 
			generator.alpha = GammaDistribution.nextGamma(2,3);
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
		treeComparison();

	}

}
