package pty.pmcmc;

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
import ma.MSAPoset;
import ma.SequenceType;
import nuts.io.CSV;
import nuts.io.IO;
import nuts.math.Sampling;
import nuts.math.StatisticsMap.DescriptiveStatisticsMap;
import nuts.util.Counter;
import pty.RandomRootedTrees;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.io.Dataset.DatasetUtils;
import pty.io.TreeEvaluator;
import pty.io.TreeEvaluator.TreeMetric;
import pty.mcmc.PhyloSampler;
import pty.mcmc.ProposalDistribution;
import pty.mcmc.UnrootedTreeState;
import pty.smc.AnnealingKernel;
import pty.smc.BackForwardKernel;
import pty.smc.LazyParticleFilter;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.LazyParticleFilter.ParticleFilterOptions;
import pty.smc.NCPriorPriorKernel;
import pty.smc.PartialCoalescentState;
import pty.smc.PartialCoalescentState4BackForwardKernel;
import pty.smc.ParticleFilter;
import pty.smc.ParticleFilterSMCSampler;
import pty.smc.ParticleFilter.ResamplingStrategy;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.ParticleFilterSMCSampler.ParticleProcessor;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.models.CTMC;
import pty.smc.test.TestBrownianModel.KernelType;
import conifer.trees.StandardNonClockPriorDensity;
import ev.ex.PhyloSamplerMain;
import ev.ex.TreeGenerators;
import ev.poi.processors.TreeDistancesProcessor;
import ev.poi.processors.TreeTopologyProcessor;
import ev.to.MrBayes;
import ev.to.NJ;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import fig.prob.Dirichlet;
import fig.prob.Gamma;
import goblin.Taxon;

public class PMCMCExperiments implements Runnable {
	@Option
	public boolean resampleRoot = false;
	@Option
	public InferenceMethod refMethod = InferenceMethod.MB;
	@Option
	public double nThousandIters = 10;
	@Option
	public ArrayList<InferenceMethod> methods = list(Arrays.asList(
			InferenceMethod.PMMHNC, InferenceMethod.SMC));
	@Option
	public ArrayList<Double> iterScalings = list(Arrays.asList(1.0));
	@Option
	public int refIterScaling = 100;
	@Option
	public int repPerDataPt = 1;
	@Option
	public double pmcmcSMCExpMix = 2.0 / 3.0;
	@Option
	public double treeRatePrior = 1.0;
	public static ev.ex.DataGenerator.DataGeneratorMain generator = new ev.ex.DataGenerator.DataGeneratorMain();
	public static PhyloSamplerMain samplerMain = new PhyloSamplerMain();
	@Option
	public static Random mainRand = new Random(3);
	@Option
	public boolean verbose = false;
	@Option
	public int nThreads = 1;
	@Option
	public int maxNUniqueParticles = Integer.MAX_VALUE;
	@Option
	public int finalMaxNUniqueParticles = Integer.MAX_VALUE;
	@Option
	public int nParticlesEachStep = 1000;
	@Option
	public double parameter_a = 1.2;
	@Option
	public double a_alpha = 1.5;
	@Option
	public double a_statFreqs = 100;
	@Option
	public double a_subsRates = 100;
	@Option
	public double a_pInv = 0.2;
	@Option
	public double pmmhPgsExpMix = 0.05;
	@Option
	public double burninPercent = 0.2; // in PMCMC algorithm
	@Option
	public int sampleTreeEveryNIter = 500;
	@Option
	public String RcommandDir = "/Users/lwa68/workspace/alpha-liangliang/Rcode/";
	// @Option public String
	// RcommandDir="/Users/l.wang/workspace/alpha-liangliang/Rcode/";
	@Option
	public boolean useDataGenerator = true;
	@Option
	public File dataFile = null;
	@Option
	public File refTree = null;
	@Option
	public String dataDirName = "output";
	@Option
	public SequenceType sequenceType = SequenceType.RNA;
	@Option
	public boolean isPMCMC4clock = true;
	@Option
	public boolean useNJinfo = false;
	@Option
	public boolean useTopologyProcessor = false;
	@Option
	public boolean saveTreesFromPMCMC = false;
	@Option
	public String nameOfAllTrees = "allTrees.trees";
	@Option
	public double csmc_trans2tranv = 2.0;
	@Option
	public double smcmcmcMix = 0.5;
	@Option
	public boolean betterStartVal = true;
	@Option
	public int nCSMC = 10;
	@Option
	public int nUCSMC = 10;
	@Option
	public boolean adaptiveTempDiff = false;
	@Option
	public double alphaSMCSampler = 0.95;

	@Option
	public double essRatioThreshold = 0.7;

	public File data = null;
	private File output = null;
	private RootedTree goldrt;

	public static void main(String[] args) {
		IO.run(args, new PMCMCExperiments(), "pcs",
				PartialCoalescentState.class, "mb", MrBayes.instance, "gen",
				generator, "mcmc", samplerMain, "prop",
				ProposalDistribution.Util._defaultProposalDistributionOptions,
				"phylo", PhyloSampler._defaultPhyloSamplerOptions, "prior",
				PhyloSampler._defaultPriorOptions, "nj", NJ.class, "nc",
				NCPriorPriorKernel.class, "priorprior", PriorPriorKernel.class,
				"annealingKernel", AnnealingKernel.class);
	}

	protected void treeComparison() {
		// data = new File( generator.output, "sim-0.msf");
		PrintWriter out = IOUtils.openOutEasy(new File(output, "results.csv"));
		out.println(CSV.header("Method", "IterScale", "Repeat", "Metric",
				"Value", "TreeName", "Time"));
		List<File> files = null;
		if (useDataGenerator)
			files = IO.ls(generator.output, "msf");
		else
			files = IO.ls(new File(Execution.getFile(dataDirName)), "msf");
		// ReportProgress.progressBlock(files.size());
		for (File f : files) {

			this.data = f;
			final String treeName = f.getName();
			LogInfo.track("Current tree:" + treeName);

			UnrootedTree goldut = (generator.useGutellData || !useDataGenerator) ? (refTree == null ? null
					: UnrootedTree.fromRooted(RootedTree.Util
							.fromNewickString(IO.f2s(refTree))))
					: UnrootedTree.fromRooted(RootedTree.Util
							.fromNewickString(IO.f2s(new File(f
									.getAbsolutePath().replaceAll("[.]msf$",
											".newick")))));

			File computedRefTrees = new File(
					Execution.getFile("computed-ref-trees"));
			// ReportProgress.progressBlock(methods.size());
			for (int i = 0; i < methods.size(); i++) {
				InferenceMethod m = methods.get(i);
				final double iterScale = iterScalings.get(i);
				LogInfo.track("Current method:" + m + " with iterScale="
						+ iterScale + " (i.e. "
						+ (iterScale * nThousandIters * 1000.0)
						+ " iterations)");

				DescriptiveStatisticsMap<String> stats = new DescriptiveStatisticsMap<String>();
				// ReportProgress.progressBlock(repPerDataPt);
				for (int j = 0; j < repPerDataPt; j++) {
					String treeNameCurrentRep = treeName;
					if (m == InferenceMethod.PMMHNC)
						treeNameCurrentRep = treeNameCurrentRep + ".Rep" + j;

					LogInfo.track("Repeat " + (j + 1) + "/" + repPerDataPt);
					// LogInfo.forceSilent = true;
					long time = System.currentTimeMillis();
					TreeDistancesProcessor processor = m.doIt(this, iterScale,
							goldut, treeNameCurrentRep);
					time = System.currentTimeMillis() - time;
					LogInfo.forceSilent = false;
					UnrootedTree inferred = processor.getConsensus();
					IO.writeToDisk(
							new File(output, "consensus_"
									+ treeName.replaceAll("[.]msf$", ".newick")),
							inferred.toNewick());
					{
						// evaluate the likelihood of the inferred tree
						Dataset dataset = DatasetUtils.fromAlignment(this.data,
								sequenceType);
						CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
						UnrootedTreeState ncs = UnrootedTreeState
								.initFastState(inferred, dataset, ctmc);
						out.println(CSV.body(m, iterScale, j, "ConsensusLogLL",
								ncs.logLikelihood(), treeName, time));
					}
					{
						// best log likelihood, when available
						double bestLogLL = processor.getBestLogLikelihood();
						out.println(CSV.body(m, iterScale, j,
								"BestSampledLogLL", bestLogLL, treeName, time));
					}
					if (goldut == null) {
						LogInfo.logsForce("Computing gold tree using "
								+ refMethod);
						goldut = refMethod.doIt(this, refIterScaling, goldut,
								treeName).getConsensus();
						computedRefTrees.mkdir();
						File current = new File(computedRefTrees,
								"computedRefTree_"
										+ f.getName().replaceAll("[.]msf", "")
										+ ".newick");
						IO.writeToDisk(current, goldut.toNewick());
					}
					for (TreeMetric tm : TreeEvaluator.coreTreeMetrics) {
						final double value = tm.score(inferred, goldut);
						stats.addValue(tm.toString(), value);

						out.println(CSV.body(m, iterScale, j, tm, value,
								treeName, time));
					}
					// if(m==InferenceMethod.PMMHNC){
					// List<Map<String,String>> pmmcRe = IO.iCSVMap(new
					// File(output, treeNameCurrentRep+".csv"));
					// Set<String> colNames=pmmcRe.get(0).keySet();
					// colNames.remove("rfDist");
					// colNames.remove("treeSize");
					// colNames.remove("Iter");
					// for(String colName:colNames)
					// {
					// DescriptiveStatisticsMap<String> statsMap=new
					// DescriptiveStatisticsMap<String>();
					// for(int k=0;k<pmmcRe.size();k++)
					// {
					// double value =
					// Double.parseDouble(pmmcRe.get(k).get(colName));
					// statsMap.addValue(colName, value);
					// }
					// out.println(CSV.body(m, iterScale, j, "pmmc_"+colName,
					// statsMap.median(colName), treeName, time));
					// }
					// }

					LogInfo.end_track();
					// ReportProgress.divisionCompleted();
				}
				LogInfo.track("Score for current block of repeats (Method=" + m
						+ ",IterScale=" + iterScale + ",TreeName=" + treeName
						+ ")");
				for (TreeMetric tm : TreeEvaluator.coreTreeMetrics)
					LogInfo.logsForce("Current " + tm + ":"
							+ stats.median(tm.toString()));
				LogInfo.end_track();
				out.flush();
				LogInfo.end_track();
				// ReportProgress.divisionCompleted();
			}
			LogInfo.end_track();

			// ReportProgress.divisionCompleted();
		}

		out.close();
	}

	public static enum InferenceMethod {
		/*
		 * MCMC {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * PhyloSampler._defaultPhyloSamplerOptions.rand = mainRand;
		 * samplerMain.alignmentInputFile = instance.data; samplerMain.st =
		 * instance.sequenceType;
		 * PhyloSampler._defaultPhyloSamplerOptions.nIteration = (int)
		 * (iterScale * instance.nThousandIters * 1000); samplerMain.run();
		 * return samplerMain.tdp; } },
		 */
		SMC {
			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilter<PartialCoalescentState> pc = new ParticleFilter<PartialCoalescentState>();
				pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				pc.nThreads = 1;
				pc.resampleLastRound = false;
				pc.rand = mainRand;

				Dataset dataset = DatasetUtils.fromAlignment(instance.data,
						instance.sequenceType);

				KernelType kernelType = KernelType.PRIOR_PRIOR;
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
				PartialCoalescentState init = PartialCoalescentState.initState(
						dataset, ctmc);
				ParticleKernel<PartialCoalescentState> ppk = kernelType.load(
						init, null);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				pc.sample(ppk, tdp);
				return tdp;
			}
		},
		MB {
			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				MrBayes mb = MrBayes.instance;
				mb.nChains = 1;
				mb.seed = mainRand.nextInt();
				mb.nMCMCIters = (int) (iterScale * instance.nThousandIters * 1000);
				if (mb.fixGTRGammaPara) {
					mb.alpha = instance.generator.alpha;
					mb.subsRates = instance.generator.subsRates;
					mb.stationaryDistribution = instance.generator.stationaryDistribution;
				}
				mb.computeSamples(MSAParser.parseMSA(instance.data),
						instance.sequenceType);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				mb.processMrBayesTrees(tdp, 1);
				// mb.cleanUpMrBayesOutput();
				return tdp;
			}
		},
		/*
		 * NCPF {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilter<PartialCoalescentState> pc = new
		 * ParticleFilter<PartialCoalescentState>(); // ParticleFilter<LazyPCS>
		 * pc = new ParticleFilter<LazyPCS>(); pc.N = (int) (iterScale *
		 * instance.nThousandIters * 1000); pc.nThreads = 8;
		 * pc.resampleLastRound = true; pc.rand = instance.mainRand; pc.verbose
		 * = instance.verbose; //
		 * pc.resamplingStrategy=ResamplingStrategy.NEVER;
		 * pc.resamplingStrategy=ResamplingStrategy.ALWAYS; Dataset dataset =
		 * DatasetUtils.fromAlignment(instance.data, instance.sequenceType);
		 * 
		 * // KernelType kernelType = KernelType.PRIOR_PRIOR;
		 * 
		 * CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
		 * PartialCoalescentState init =
		 * PartialCoalescentState.initFastState(dataset, ctmc, false);
		 * ParticleKernel<PartialCoalescentState> ppk = new
		 * NCPriorPriorKernel(init); TreeDistancesProcessor tdp = new
		 * TreeDistancesProcessor();
		 * 
		 * pc.sample(ppk,tdp);
		 * System.out.println("initial loglikelihood: "+init.
		 * logLikelihood()+" Marginal loglikelihood: "
		 * +(pc.estimateNormalizer()+init.logLikelihood())); return tdp; } },
		 * LAZYNCPF {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilter<LazyPCS> pc = new ParticleFilter<LazyPCS>(); pc.N =
		 * (int) (iterScale * instance.nThousandIters * 1000); pc.nThreads = 1;
		 * pc.resampleLastRound = true; pc.rand = instance.mainRand; Dataset
		 * dataset = DatasetUtils.fromAlignment(instance.data,
		 * instance.sequenceType); CTMC ctmc =
		 * CTMC.SimpleCTMC.dnaCTMC(dataset.nSites()); PartialCoalescentState
		 * init = PartialCoalescentState.initFastState(dataset, ctmc, false);
		 * LazyNonclockPriorPrior ppk = new LazyNonclockPriorPrior(init);
		 * TreeDistancesProcessor tdp = new TreeDistancesProcessor();
		 * pc.sample(ppk,tdp);
		 * System.out.println("Marginal loglikelihood: "+pc.estimateNormalizer
		 * ()); return tdp; } },
		 */
		SMC4K2P {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {

				// if (instance.usePriorPost)
				// throw new RuntimeException();
				//
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

				Dataset dataset = DatasetUtils.fromAlignment(instance.data,
						instance.sequenceType);
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
				PartialCoalescentState init = PartialCoalescentState
						.initFastState(dataset, ctmc);
				LazyParticleKernel pk2 = new PriorPriorKernel(init);
				LazyParticleFilter<PartialCoalescentState> lpf = new LazyParticleFilter<PartialCoalescentState>(
						pk2, options);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();

				if (instance.useTopologyProcessor) {
					TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
					final double zHat = lpf.sample(tdp, trTopo);

					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: "
							+ urtCounter.keySet().size());
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
		SMCNonClock4K2P {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {

				// if (instance.usePriorPost)
				// throw new RuntimeException();

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

				Dataset dataset = DatasetUtils.fromAlignment(instance.data,
						instance.sequenceType);
				// CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),
						instance.csmc_trans2tranv);
				PartialCoalescentState init = PartialCoalescentState
						.initFastState(instance.resampleRoot, dataset, ctmc,
								false);
				LazyParticleKernel pk2 = new NCPriorPriorKernel(init);
				LazyParticleFilter<PartialCoalescentState> lpf = new LazyParticleFilter<PartialCoalescentState>(
						pk2, options);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				final double zHat = lpf.sample(tdp);
				LogInfo.logsForce("Norm:" + zHat);
				return tdp;
			}
		},
		MIXSMCNCMB {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilterOptions options = new ParticleFilterOptions();
				options.nParticles = (int) (iterScale * instance.nThousandIters
						* 1000 * instance.smcmcmcMix);
				// options.nThreads = instance.nThreads;
				options.nThreads = 1;
				options.resampleLastRound = true;
				options.parallelizeFinalParticleProcessing = true;
				options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
				options.maxNUniqueParticles = instance.maxNUniqueParticles;
				options.rand = instance.mainRand;
				options.verbose = instance.verbose;
				Dataset dataset = DatasetUtils.fromAlignment(instance.data,
						instance.sequenceType);
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),
						instance.csmc_trans2tranv);
				PartialCoalescentState init = PartialCoalescentState
						.initFastState(instance.resampleRoot, dataset, ctmc,
								false);
				LazyParticleKernel pk2 = new NCPriorPriorKernel(init);
				LazyParticleFilter<PartialCoalescentState> lpf = new LazyParticleFilter<PartialCoalescentState>(
						pk2, options);
				StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();
				// TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				final double zHat = lpf.sample(pro);
				LogInfo.logsForce("Norm:" + zHat);
				PartialCoalescentState sampled = pro.sample(options.rand);
				RootedTree rtSample = sampled.getFullCoalescentState();
				String stringOfTree = RootedTree.Util.toNewick(rtSample);
				MrBayes mb = MrBayes.instance;
				mb.nChains = 1;
				mb.seed = mainRand.nextInt();
				mb.nMCMCIters = (int) (iterScale * instance.nThousandIters
						* 1000 * (1 - instance.smcmcmcMix));
				mb.setStartTree(stringOfTree);
				if (mb.fixGTRGammaPara) {
					mb.alpha = instance.generator.alpha;
					mb.subsRates = instance.generator.subsRates;
					mb.stationaryDistribution = instance.generator.stationaryDistribution;
				}
				mb.computeSamples(MSAParser.parseMSA(instance.data),
						instance.sequenceType);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				mb.processMrBayesTrees(tdp, 1);
				return tdp;
			}
		},
		/*
		 * // For nonclock tree PMMHSMC2NC {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilterOptions options = new ParticleFilterOptions();
		 * options.nParticles = instance.nParticlesEachStep; //
		 * generator.nTaxa*(generator.nTaxa-1)*20; //(int) (iterScale *
		 * instance.nThousandIters * 1000); options.nThreads =
		 * instance.nThreads; options.resampleLastRound = true;
		 * options.parallelizeFinalParticleProcessing = true;
		 * options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
		 * options.maxNUniqueParticles = instance.maxNUniqueParticles;
		 * options.rand = instance.mainRand; options.verbose = instance.verbose;
		 * // options.maxNGrow = 0;
		 * 
		 * Dataset dataset = DatasetUtils.fromAlignment(instance.data,
		 * instance.sequenceType); CTMC ctmc =
		 * CTMC.SimpleCTMC.dnaCTMC(dataset.nSites()); PartialCoalescentState
		 * init = PartialCoalescentState.initFastState(dataset, ctmc,false);
		 * LazyParticleKernel pk2 = new NCPriorPriorKernel(init);
		 * LazyParticleFilter<PartialCoalescentState> lpf = new
		 * LazyParticleFilter<PartialCoalescentState>(pk2, options);
		 * TreeDistancesProcessor tdp = new TreeDistancesProcessor(); double
		 * initTrans2tranv=
		 * Math.exp(Math.log(2.0)+instance.mainRand.nextGaussian()*0.1);
		 * PMMHSMC2NC pmmh = new PMMHSMC2NC(dataset, lpf, tdp, initTrans2tranv);
		 * pmmh.a=instance.parameter_a; // final double requested = (iterScale *
		 * instance.nThousandIters * 1000); // pc.N = 1+ (int)
		 * Math.pow(requested, instance.pmcmcSMCExpMix); // final int nMCMC = 1
		 * + (int) Math.pow(requested, 1.0 - instance.pmcmcSMCExpMix); final int
		 * nMCMC=(int) iterScale; System.out.println("# particles: "+
		 * options.nParticles +"; nMCMC:"+ nMCMC); int i=0; // for (int i = 0; i
		 * < nMCMC; i++) while(!pmmh.stop && i<nMCMC) { i++;
		 * pmmh.next(instance.mainRand); } return tdp; } }, // For clock tree
		 * PMMHSMC2 {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilterOptions options = new ParticleFilterOptions();
		 * options.nParticles = instance.nParticlesEachStep; //
		 * generator.nTaxa*(generator.nTaxa-1)*20; //(int) (iterScale *
		 * instance.nThousandIters * 1000); options.nThreads =
		 * instance.nThreads; options.resampleLastRound = true;
		 * options.parallelizeFinalParticleProcessing = true;
		 * options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
		 * options.maxNUniqueParticles = instance.maxNUniqueParticles;
		 * options.rand = instance.mainRand; options.verbose = instance.verbose;
		 * // options.maxNGrow = 0;
		 * 
		 * Dataset dataset = DatasetUtils.fromAlignment(instance.data,
		 * instance.sequenceType); CTMC ctmc =
		 * CTMC.SimpleCTMC.dnaCTMC(dataset.nSites()); PartialCoalescentState
		 * init = PartialCoalescentState.initFastState(dataset, ctmc);
		 * LazyParticleKernel pk2 = new PriorPriorKernel(init);
		 * LazyParticleFilter<PartialCoalescentState> lpf = new
		 * LazyParticleFilter<PartialCoalescentState>(pk2, options);
		 * TreeDistancesProcessor tdp = new TreeDistancesProcessor(); double
		 * initTrans2tranv=
		 * Math.exp(Math.log(2.0)+instance.mainRand.nextGaussian()*0.1);
		 * PMMHSMC2 pmmh = new PMMHSMC2(dataset, options, tdp, initTrans2tranv);
		 * pmmh.a=instance.parameter_a; // final double requested = (iterScale *
		 * instance.nThousandIters * 1000); // pc.N = 1+ (int)
		 * Math.pow(requested, instance.pmcmcSMCExpMix); // final int nMCMC = 1
		 * + (int) Math.pow(requested, 1.0 - instance.pmcmcSMCExpMix); final int
		 * nMCMC=(int) iterScale; System.out.println("# particles: "+
		 * options.nParticles +"; nMCMC:"+ nMCMC); int i=0; // for (int i = 0; i
		 * < nMCMC; i++) while(!pmmh.stop && i<nMCMC) { i++;
		 * pmmh.next(instance.mainRand); } return tdp; } },
		 */

		PMMHNC {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilter<PartialCoalescentState> pc = new ParticleFilter<PartialCoalescentState>();
				// pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				pc.N = instance.nParticlesEachStep; // generator.nTaxa*(generator.nTaxa-1)*100;
				pc.nThreads = instance.nThreads;
				pc.resampleLastRound = false;
				pc.resamplingStrategy = ResamplingStrategy.ALWAYS;
				pc.rand = instance.mainRand;
				Dataset dataset = DatasetUtils.fromAlignment(instance.data,
						instance.sequenceType);
				double initTrans2tranv = Math.exp(Math.log(2.0)
						+ instance.mainRand.nextGaussian() * 0.1);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				PMMHNC pmmhnc = new PMMHNC(dataset, pc, tdp, initTrans2tranv);
				// pmmhnc.goldut=goldut;
				// pmmhnc.out=treeName;
				pmmhnc.a = instance.parameter_a;
				// final double requested = (iterScale * instance.nThousandIters
				// * 1000);
				// pc.N = 1+ (int) Math.pow(requested, instance.pmcmcSMCExpMix);
				// final int nMCMC = 1 + (int) Math.pow(requested, 1.0 -
				// instance.pmcmcSMCExpMix);
				final int nMCMC = (int) (iterScale * instance.nThousandIters * 1000);
				System.out.println("# particles: " + pc.N + "; nMCMC:" + nMCMC);
				int i = 0;
				// for (int i = 0; i < nMCMC; i++)
				while (!pmmhnc.stop && i < nMCMC) {
					i++;
					pmmhnc.next(instance.mainRand);
				}
				return tdp;
			}
		},
		SMC4GTRIGammaBF {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
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
				double[] subsRates = instance.generator.subsRates; // six
				// parameters
				// of
				// substitutions:rAC,rAG,rAT,rCG,rGT,rCT
				double[] statFreqs = instance.generator.stationaryDistribution; // stationary
				// state
				// frequencies.
				// pi_A,
				// pi_C,
				// pi_G,
				// pi_T
				double alpha = instance.generator.alpha; // shape parameter in
				// the Gamma
				// distribution
				double pInv = instance.generator.pInv; // the proportion of
				// invariant sites
				int nCategories = 4;
				int dataRepeatN = nCategories;
				if (pInv > 0)
					dataRepeatN = nCategories + 1;
				MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
				Dataset dataset = DatasetUtils.fromAlignment(align,
						instance.sequenceType, dataRepeatN);
				CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
						dataset.nSites(), alpha, nCategories, pInv);
				PartialCoalescentState init0 = PartialCoalescentState
						.initFastState(dataset, ctmc, true);
				PartialCoalescentState4BackForwardKernel init = new PartialCoalescentState4BackForwardKernel(
						init0, null, 0);
				LazyParticleKernel pk2 = new BackForwardKernel(init);
				LazyParticleFilter<PartialCoalescentState4BackForwardKernel> lpf = new LazyParticleFilter<PartialCoalescentState4BackForwardKernel>(
						pk2, options);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				if (instance.useTopologyProcessor) {
					TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
					final double zHat = lpf.sample(tdp, trTopo);
					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: "
							+ urtCounter.keySet().size());
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
		SMC4GTRIGamma {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
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
				double[] subsRates = instance.generator.subsRates; // six
																	// parameters
																	// of
																	// substitutions:rAC,rAG,rAT,rCG,rGT,rCT
				double[] statFreqs = instance.generator.stationaryDistribution; // stationary
																				// state
																				// frequencies.
																				// pi_A,
																				// pi_C,
																				// pi_G,
																				// pi_T
				double alpha = instance.generator.alpha; // shape parameter in
															// the Gamma
															// distribution
				double pInv = instance.generator.pInv; // the proportion of
														// invariant sites
				int nCategories = 4;
				int dataRepeatN = nCategories;
				if (pInv > 0)
					dataRepeatN = nCategories + 1;
				MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
				Dataset dataset = DatasetUtils.fromAlignment(align,
						instance.sequenceType, dataRepeatN);
				CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
						dataset.nSites(), alpha, nCategories, pInv);
				PartialCoalescentState init = PartialCoalescentState
						.initFastState(dataset, ctmc, true);
				LazyParticleKernel pk2 = new PriorPriorKernel(init);
				LazyParticleFilter<PartialCoalescentState> lpf = new LazyParticleFilter<PartialCoalescentState>(
						pk2, options);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				if (instance.useTopologyProcessor) {
					TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
					final double zHat = lpf.sample(tdp, trTopo);
					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: "
							+ urtCounter.keySet().size());
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
		/*
		 * SMC4GTRIGammaSigmaTheta {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilterOptions options = new ParticleFilterOptions(); //
		 * ParticleFilter<PartialCoalescentState> pc = new
		 * ParticleFilter<PartialCoalescentState>(); options.nParticles = (int)
		 * (iterScale * instance.nThousandIters * 1000); options.nThreads =
		 * instance.nThreads; options.resampleLastRound = true;
		 * options.parallelizeFinalParticleProcessing = true;
		 * options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
		 * options.maxNUniqueParticles = instance.maxNUniqueParticles;
		 * //options.rand = instance.mainRand; options.verbose =
		 * instance.verbose; // options.maxNGrow = 0; // double[] subsRates=new
		 * double[]{0.26,0.18,0.17,0.15,0.11,0.13}; // six parameters of
		 * substitutions:rAC,rAG,rAT,rCG,rGT,rCT // double[] statFreqs=new
		 * double[]{0.3,0.2,0.2,0.3}; // stationary state frequencies. pi_A,
		 * pi_C, pi_G, pi_T
		 * 
		 * double[] subsRates=instance.generator.subsRates; double[]
		 * statFreqs=instance.generator.stationaryDistribution; double
		 * alpha=instance.generator.alpha; // shape parameter in the Gamma
		 * distribution double pInv=generator.alpha; // the proportion of
		 * invariant sites int nCategories=4; int dataRepeatN=nCategories;
		 * if(pInv>0)dataRepeatN=nCategories+1; MSAPoset align =
		 * MSAPoset.parseAlnOrMsfFormats(instance.data); Dataset dataset =
		 * DatasetUtils.fromAlignment(align, instance.sequenceType,dataRepeatN);
		 * CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
		 * dataset.nSites(), alpha, nCategories, pInv); PartialCoalescentState
		 * init = PartialCoalescentState.initFastState(dataset, ctmc, true);
		 * TreeDistancesProcessor tdp=null; int M=5; double [] logLikes=new
		 * double[M]; double mu=0; //mean for(int i=0;i<M;i++) {
		 * System.out.println(i); options.rand = new
		 * Random(instance.mainRand.nextLong()); LazyParticleKernel pk2 = new
		 * PriorPriorKernel(init); LazyParticleFilter<PartialCoalescentState>
		 * lpf = new LazyParticleFilter<PartialCoalescentState>(pk2, options);
		 * tdp = new TreeDistancesProcessor(); final double zHat =
		 * lpf.sample(tdp); // LogInfo.logsForce("Norm:" + zHat);
		 * logLikes[i]=zHat; mu+=logLikes[i]; } mu=mu/M; //
		 * LogInfo.logsForce("Norm:" + zHat); double sigma2=0; for(int
		 * i=0;i<M;i++){ sigma2+=Math.pow(logLikes[i]-mu, 2);
		 * LogInfo.logsForce(logLikes[i]+","); } sigma2=sigma2/M;
		 * LogInfo.logsForce("Sample variance is "+sigma2);
		 * LogInfo.logsForce("Sample SD is "+Math.sqrt(sigma2)); return tdp; }
		 * },
		 */
		SMCNonClock4GTRIGamma {
			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
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
				double[] subsRates = instance.generator.subsRates;
				double[] statFreqs = instance.generator.stationaryDistribution;

				double alpha = instance.generator.alpha; // shape parameter in
															// the Gamma
															// distribution
				double pInv = instance.generator.pInv; // the proportion of
														// invariant sites
				int nCategories = 4;
				int dataRepeatN = nCategories;
				if (pInv > 0)
					dataRepeatN = nCategories + 1;
				MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
				Dataset dataset = DatasetUtils.fromAlignment(align,
						instance.sequenceType, dataRepeatN);
				CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
						dataset.nSites(), alpha, nCategories, pInv);
				PartialCoalescentState init = PartialCoalescentState
						.initFastState(dataset, ctmc, false);
				LazyParticleKernel pk2 = new NCPriorPriorKernel(init);
				LazyParticleFilter<PartialCoalescentState> lpf = new LazyParticleFilter<PartialCoalescentState>(
						pk2, options);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				if (instance.useTopologyProcessor) {
					TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
					final double zHat = lpf.sample(tdp, trTopo);

					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: "
							+ urtCounter.keySet().size());
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
		// SMCNonClock4GTRIGammaSamples {
		// @Override
		// public TreeDistancesProcessor doIt(PMCMCExperiments instance, double
		// iterScale, UnrootedTree goldut, String treeName)
		// {
		// ParticleFilterOptions options = new ParticleFilterOptions();
		// // ParticleFilter<PartialCoalescentState> pc = new
		// ParticleFilter<PartialCoalescentState>();
		// options.nParticles = (int) (iterScale * instance.nThousandIters *
		// 1000);
		// options.nThreads = instance.nThreads;
		// options.resampleLastRound = true;
		// options.parallelizeFinalParticleProcessing = true;
		// options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
		// options.maxNUniqueParticles = instance.maxNUniqueParticles;
		// options.rand = instance.mainRand;
		// options.verbose = instance.verbose;
		// // options.maxNGrow = 0;
		// File output=new File(Execution.getFile("results"));
		// String nameOfAllTrees="allTrees.trees";
		// double[] subsRates=instance.generator.subsRates;
		// double[] statFreqs=instance.generator.stationaryDistribution;
		//
		// double alpha=instance.generator.alpha; // shape parameter in the
		// Gamma distribution
		// double pInv=instance.generator.pInv; // the proportion of invariant
		// sites
		// int nCategories=4;
		// int dataRepeatN=nCategories;
		// if(pInv>0)dataRepeatN=nCategories+1;
		// MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
		// Dataset dataset = DatasetUtils.fromAlignment(align,
		// instance.sequenceType,dataRepeatN);
		// CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
		// dataset.nSites(), alpha, nCategories, pInv);
		// PartialCoalescentState init =
		// PartialCoalescentState.initFastState(dataset, ctmc, false);
		// LazyParticleKernel pk2 = new NCPriorPriorKernel(init);
		// LazyParticleFilter<PartialCoalescentState> lpf = new
		// LazyParticleFilter<PartialCoalescentState>(pk2, options);
		//
		// StoreProcessor<PartialCoalescentState> pro = new
		// StoreProcessor<PartialCoalescentState>();
		// TreeDistancesProcessor tdp = new TreeDistancesProcessor();
		// final double zHat = lpf.sample(pro,tdp);
		// // LogInfo.logsForce("Norm:" + zHat);
		// //int nSamples=Math.max(options.nParticles, 10000);
		// int nSamples=options.nParticles;
		// for(int i=0;i<nSamples;i++)
		// {
		// PartialCoalescentState sampled = pro.sample(instance.mainRand);
		// //TODO: there is a problem here when the nThread >1
		// RootedTree currentSample = sampled.getFullCoalescentState();
		//
		// String stringOfTree=RootedTree.Util.toNewick(currentSample);
		// String cmdStr="echo -n 'TREE tree_" + i +"=' >>"+nameOfAllTrees;
		// IO.call("bash -s",cmdStr,output);
		// cmdStr="echo '" +stringOfTree +
		// "' | sed 's/internal_[0-9]*_[0-9]*//g'"+ " | sed 's/leaf_//g'"+
		// " >> "+nameOfAllTrees;
		// // LogInfo.logs(cmdStr);
		// // String msg0=
		// IO.call("bash -s",cmdStr,output);
		// }
		// return tdp;
		// }
		// },
		/*
		 * SMCNonClock4GTRIGammaSigmaTheta {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilterOptions options = new ParticleFilterOptions(); //
		 * ParticleFilter<PartialCoalescentState> pc = new
		 * ParticleFilter<PartialCoalescentState>(); options.nParticles = (int)
		 * (iterScale * instance.nThousandIters * 1000); options.nThreads =
		 * instance.nThreads; options.resampleLastRound = true;
		 * options.parallelizeFinalParticleProcessing = true;
		 * options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
		 * options.maxNUniqueParticles = instance.maxNUniqueParticles;
		 * //options.rand = instance.mainRand; options.verbose =
		 * instance.verbose; // options.maxNGrow = 0; double[]
		 * subsRates=instance.generator.subsRates; double[]
		 * statFreqs=instance.generator.stationaryDistribution;
		 * 
		 * double alpha=instance.generator.alpha; // shape parameter in the
		 * Gamma distribution double pInv=instance.generator.pInv; // the
		 * proportion of invariant sites int nCategories=4; int
		 * dataRepeatN=nCategories; if(pInv>0)dataRepeatN=nCategories+1;
		 * MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
		 * Dataset dataset = DatasetUtils.fromAlignment(align,
		 * instance.sequenceType,dataRepeatN); CTMC ctmc = new
		 * CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4, dataset.nSites(), alpha,
		 * nCategories, pInv); PartialCoalescentState init =
		 * PartialCoalescentState.initFastState(dataset, ctmc, false);
		 * LazyParticleKernel pk2 = new NCPriorPriorKernel(init);
		 * TreeDistancesProcessor tdp=null; int M=5; double [] logLikes=new
		 * double[M]; double mu=0; //mean for(int i=0;i<M;i++) {
		 * System.out.println(i); options.rand = new
		 * Random(instance.mainRand.nextLong()); //new Random(instance.mainRand)
		 * LazyParticleFilter<PartialCoalescentState> lpf = new
		 * LazyParticleFilter<PartialCoalescentState>(pk2, options); tdp = new
		 * TreeDistancesProcessor(); final double zHat = lpf.sample(tdp);
		 * logLikes[i]=zHat; mu+=logLikes[i]; } mu=mu/M; //
		 * LogInfo.logsForce("Norm:" + zHat); double sigma2=0; for(int
		 * i=0;i<M;i++){ sigma2+=Math.pow(logLikes[i]-mu, 2);
		 * LogInfo.logsForce(logLikes[i]+","); } sigma2=sigma2/M;
		 * LogInfo.logsForce("Sample variance is "+sigma2);
		 * LogInfo.logsForce("Sample SD is "+Math.sqrt(sigma2)); return tdp; }
		 * }, SMC_NJ_NonClock4GTRIGamma {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilterOptions options = new ParticleFilterOptions(); //
		 * ParticleFilter<PartialCoalescentState> pc = new
		 * ParticleFilter<PartialCoalescentState>(); options.nParticles = (int)
		 * (iterScale * instance.nThousandIters * 1000); options.nThreads =
		 * instance.nThreads; options.resampleLastRound = true;
		 * options.parallelizeFinalParticleProcessing = true;
		 * options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
		 * options.maxNUniqueParticles = instance.maxNUniqueParticles;
		 * options.rand = instance.mainRand; options.verbose = instance.verbose;
		 * // options.maxNGrow = 0; double[]
		 * subsRates=instance.generator.subsRates; double[]
		 * statFreqs=instance.generator.stationaryDistribution; double
		 * alpha=instance.generator.alpha; // shape parameter in the Gamma
		 * distribution double pInv=instance.generator.alpha; // the proportion
		 * of invariant sites // int nCategories=4; // int
		 * dataRepeatN=nCategories; // if(pInv>0)dataRepeatN=nCategories+1; //
		 * MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data); //
		 * Dataset dataset = DatasetUtils.fromAlignment(align,
		 * instance.sequenceType,dataRepeatN); // CTMC ctmc = new
		 * CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4, dataset.nSites(), alpha,
		 * nCategories, pInv); //
		 * 
		 * Dataset dataset = DatasetUtils.fromAlignment(instance.data,
		 * instance.sequenceType); CTMC ctmc =
		 * CTMC.SimpleCTMC.dnaCTMC(dataset.nSites()); PartialCoalescentState
		 * init = PartialCoalescentState.initFastState(instance.resampleRoot,
		 * dataset, ctmc, false);
		 * 
		 * 
		 * //PartialCoalescentState init =
		 * PartialCoalescentState.initFastState(dataset, ctmc, false); NJPState2
		 * initNJP = NJPState2.init(init); // LazyParticleKernel pk2 = new
		 * NCPriorPriorKernel(init);
		 * 
		 * // LazyParticleKernel pk2 =new LazyNJStateNonclockKernel(initNJP);
		 * 
		 * LazyParticleKernel pk2 =new NJStateNonclockKernel(initNJP);
		 * LazyParticleFilter<PartialCoalescentState> lpf = new
		 * LazyParticleFilter<PartialCoalescentState>(pk2, options);
		 * TreeDistancesProcessor tdp = new TreeDistancesProcessor(); final
		 * double zHat = lpf.sample(tdp); // LogInfo.logsForce("Norm:" + zHat);
		 * return tdp; } },
		 */
		PGS4GTRIGammaBF {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilterOptions options = new ParticleFilterOptions();
				options.nParticles = instance.nParticlesEachStep; // generator.nTaxa*(generator.nTaxa-1)*20;
																	// //(int)
																	// (iterScale
																	// *
																	// instance.nThousandIters
																	// * 1000);
				// options.nThreads = instance.nThreads;
				options.nThreads = 1; // TODO: solve the problems of using
										// multiple threads in pmmh.
				options.resampleLastRound = true;
				options.parallelizeFinalParticleProcessing = true;
				options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
				options.maxNUniqueParticles = instance.maxNUniqueParticles;
				options.rand = instance.mainRand;
				options.verbose = instance.verbose;
				// six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
				double[] subsRates = Dirichlet.sample(instance.mainRand,
						new double[] { 10, 10, 10, 10, 10, 10 });
				// double[] subsRates=new
				// double[]{0.26,0.18,0.17,0.15,0.11,0.13};
				// stationary state frequencies. pi_A, pi_C, pi_G, pi_T
				double[] statFreqs = Dirichlet.sample(instance.mainRand,
						new double[] { 10, 10, 10, 10 });
				// double[] statFreqs=new double[]{0.3,0.2,0.2,0.3};
				double alpha = Sampling.nextDouble(instance.mainRand, 0.1, 0.9); // shape
																					// parameter
																					// in
																					// the
																					// Gamma
																					// distribution
				if (instance.betterStartVal) {
					double[] subsRatesTmp = new double[subsRates.length];
					for (int i = 0; i < subsRates.length; i++)
						subsRatesTmp[i] = instance.generator.subsRates[i] * 10;
					subsRates = Dirichlet.sample(instance.mainRand,
							subsRatesTmp);
					double[] subStatFreqsTmp = new double[statFreqs.length];
					for (int i = 0; i < statFreqs.length; i++)
						subStatFreqsTmp[i] = instance.generator.stationaryDistribution[i] * 10;
					statFreqs = Dirichlet.sample(instance.mainRand,
							subStatFreqsTmp);
					alpha = Sampling.nextDouble(instance.mainRand,
							Math.max(instance.generator.alpha - 0.1, 0),
							Math.max(instance.generator.alpha + 0.1, 1.0));
				}

				double pInv = 0; // the proportion of invariant sites
				int nCategories = 4;
				int dataRepeatN = nCategories;
				if (pInv > 0)
					dataRepeatN = nCategories + 1;
				MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
				Dataset dataset = DatasetUtils.fromAlignment(align,
						instance.sequenceType, dataRepeatN);
				// CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
				// dataset.nSites(), alpha, nCategories, pInv);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
				final int nMCMC = (int) iterScale;
				final int nPMMH = 0;
				final int nPGS = nMCMC - nPMMH;
				String resultFolder = Execution.getFile("results");
				File output = new File(resultFolder);
				LogInfo.logsForce(" # particles: " + options.nParticles
						+ "; nMCMC:" + nMCMC);
				RootedTree initTree = null;
				if (nPMMH == 0)
					initTree = RandomRootedTrees.sampleCoalescent(
							instance.mainRand, align.nTaxa(), 10);
				options.nThreads = instance.nThreads;
				ParticleGibbs4GTRIGamma pg = new ParticleGibbs4GTRIGamma(
						dataset, options, tdp, instance.useTopologyProcessor,
						trTopo, initTree, subsRates, statFreqs, alpha, pInv,
						instance.a_alpha, instance.a_pInv,
						instance.a_statFreqs, instance.a_subsRates,
						nCategories, false, instance.isPMCMC4clock,
						instance.sampleTreeEveryNIter);
				pg.setSaveTreesFromPMCMC(instance.saveTreesFromPMCMC);
				pg.setNameOfAllTrees(instance.nameOfAllTrees);
				int i = 0;
				final int nPGSburnin = (int) (nPGS * instance.burninPercent);
				while (i < nPGS) {
					i++;
					System.out.println(i);
					if (i > nPGSburnin)
						pg.setProcessTree(true);
					pg.next(instance.mainRand);
				}
				if (instance.saveTreesFromPMCMC) {
					IO.call("bash -s", "echo 'END;' >> "
							+ instance.nameOfAllTrees, output);
					String RcmdStr = instance.RcommandDir
							+ "thetaTracePlot.R  \"resultFolder='"
							+ resultFolder + "'\"";
					String msg0 = IO.call("bash -s", RcmdStr, output);
					LogInfo.logs(msg0);
				}

				if (instance.useTopologyProcessor) {
					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: "
							+ urtCounter.keySet().size());
					for (UnrootedTree urt : urtCounter.keySet()) {
						LogInfo.logsForce(urt);
						LogInfo.logsForce(urtCounter.getCount(urt));
					}
				}
				return tdp;
			}
		},
		PMMH_PGS4GTRIGamma {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilterOptions options = new ParticleFilterOptions();
				options.nParticles = instance.nParticlesEachStep; // generator.nTaxa*(generator.nTaxa-1)*20;
				// //(int)
				// (iterScale
				// *
				// instance.nThousandIters
				// * 1000);
				// options.nThreads = instance.nThreads;
				options.nThreads = 1; // TODO: solve the problems of using
				// multiple threads in pmmh.
				options.resampleLastRound = true;
				options.parallelizeFinalParticleProcessing = true;
				options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
				options.maxNUniqueParticles = instance.maxNUniqueParticles;
				options.rand = instance.mainRand;
				options.verbose = instance.verbose;
				// six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
				double[] subsRates = Dirichlet.sample(instance.mainRand,
						new double[] { 10, 10, 10, 10, 10, 10 });
				// double[] subsRates=new
				// double[]{0.26,0.18,0.17,0.15,0.11,0.13};
				// stationary state frequencies. pi_A, pi_C, pi_G, pi_T
				double[] statFreqs = Dirichlet.sample(instance.mainRand,
						new double[] { 10, 10, 10, 10 });
				// double[] statFreqs=new double[]{0.3,0.2,0.2,0.3};
				double alpha = Sampling.nextDouble(instance.mainRand, 0.1, 0.9); // shape
				// parameter
				// in
				// the
				// Gamma
				// distribution
				if (instance.betterStartVal) {
					double[] subsRatesTmp = new double[subsRates.length];
					for (int i = 0; i < subsRates.length; i++)
						subsRatesTmp[i] = instance.generator.subsRates[i] * 10;
					subsRates = Dirichlet.sample(instance.mainRand,
							subsRatesTmp);
					double[] subStatFreqsTmp = new double[statFreqs.length];
					for (int i = 0; i < statFreqs.length; i++)
						subStatFreqsTmp[i] = instance.generator.stationaryDistribution[i] * 10;
					statFreqs = Dirichlet.sample(instance.mainRand,
							subStatFreqsTmp);
					alpha = Sampling.nextDouble(instance.mainRand,
							Math.max(instance.generator.alpha - 0.1, 0),
							Math.max(instance.generator.alpha + 0.1, 1.0));
				}

				double pInv = 0; // the proportion of invariant sites
				int nCategories = 4;
				int dataRepeatN = nCategories;
				if (pInv > 0)
					dataRepeatN = nCategories + 1;
				MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
				Dataset dataset = DatasetUtils.fromAlignment(align,
						instance.sequenceType, dataRepeatN);
				// CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
				// dataset.nSites(), alpha, nCategories, pInv);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
				final int nMCMC = (int) iterScale;
				final int nPMMH = (int) (nMCMC * instance.pmmhPgsExpMix);
				final int nPGS = nMCMC - nPMMH;
				String resultFolder = Execution.getFile("results");
				File output = new File(resultFolder);
				PMMH4GTRIGamma pmmh = new PMMH4GTRIGamma(dataset, options, tdp,
						instance.useTopologyProcessor, trTopo, subsRates,
						statFreqs, alpha, pInv, instance.a_alpha,
						instance.a_pInv, instance.a_statFreqs,
						instance.a_subsRates, nCategories, false,
						instance.isPMCMC4clock, instance.useNJinfo);
				pmmh.setSaveTreesFromPMCMC(instance.saveTreesFromPMCMC);
				if (instance.saveTreesFromPMCMC) {
					String cmdStr = "echo '#NEXUS' > "
							+ instance.nameOfAllTrees;
					String msg0 = IO.call("bash -s", cmdStr, output);
					IO.call("bash -s", "echo 'BEGIN TREES;' >> "
							+ instance.nameOfAllTrees, output);
					IO.call("bash -s", "echo 'TRANSLATE' >> "
							+ instance.nameOfAllTrees, output);
					for (int j = 0; j < dataset.observations().size(); j++)
						IO.call("bash -s", "echo '" + j + " leaf_" + j
								+ ",' >> " + instance.nameOfAllTrees, output);
					IO.call("bash -s",
							"echo ';' >> " + instance.nameOfAllTrees, output);
					pmmh.setNameOfAllTrees(instance.nameOfAllTrees);
				}
				LogInfo.logsForce(" # particles: " + options.nParticles
						+ "; nMCMC:" + nMCMC);
				final int nPMMHburnin = (int) (nPMMH * instance.burninPercent);
				int i = 0;
				while (i < nPMMH) {
					i++;
					// System.out.println(i);
					if (i > nPMMHburnin)
						pmmh.setProcessTree(true);
					pmmh.next(instance.mainRand);
				}
				RootedTree initTree = pmmh.getRootedTree();
				if (nPMMH == 0)
					initTree = RandomRootedTrees.sampleCoalescent(
							instance.mainRand, align.nTaxa(), 10);
				options.nThreads = instance.nThreads;
				ParticleGibbs4GTRIGamma pg = new ParticleGibbs4GTRIGamma(
						dataset, options, tdp, instance.useTopologyProcessor,
						trTopo, initTree, pmmh.getSubsRates(),
						pmmh.getStateFreqs(), pmmh.getAlpha(), pmmh.getpInv(),
						instance.a_alpha, instance.a_pInv,
						instance.a_statFreqs, instance.a_subsRates,
						nCategories, false, instance.isPMCMC4clock,
						instance.sampleTreeEveryNIter);
				pg.setSaveTreesFromPMCMC(instance.saveTreesFromPMCMC);
				pg.setNameOfAllTrees(instance.nameOfAllTrees);
				i = 0;
				final int nPGSburnin = (int) (nPGS * instance.burninPercent);
				while (i < nPGS) {
					i++;
					System.out.println(i);
					if (i > nPGSburnin)
						pg.setProcessTree(true);
					pg.next(instance.mainRand);
				}
				if (instance.saveTreesFromPMCMC) {
					IO.call("bash -s", "echo 'END;' >> "
							+ instance.nameOfAllTrees, output);
					String RcmdStr = instance.RcommandDir
							+ "thetaTracePlot.R  \"resultFolder='"
							+ resultFolder + "'\"";
					String msg0 = IO.call("bash -s", RcmdStr, output);
					LogInfo.logs(msg0);
				}

				if (instance.useTopologyProcessor) {
					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: "
							+ urtCounter.keySet().size());
					for (UnrootedTree urt : urtCounter.keySet()) {
						LogInfo.logsForce(urt);
						LogInfo.logsForce(urtCounter.getCount(urt));
					}
				}
				return tdp;
			}
		},
		PMMH_iPGS4GTRIGamma {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilterOptions options = new ParticleFilterOptions();
				options.nParticles = instance.nParticlesEachStep; // generator.nTaxa*(generator.nTaxa-1)*20;
				// //(int)
				// (iterScale
				// *
				// instance.nThousandIters
				// * 1000);
				// options.nThreads = instance.nThreads;
				options.nThreads = 1; // TODO: solve the problems of using
				// multiple threads in pmmh.
				options.resampleLastRound = true;
				options.parallelizeFinalParticleProcessing = true;
				options.finalMaxNUniqueParticles = instance.finalMaxNUniqueParticles;
				options.maxNUniqueParticles = instance.maxNUniqueParticles;
				options.rand = instance.mainRand;
				options.verbose = instance.verbose;
				// six parameters of substitutions:rAC,rAG,rAT,rCG,rGT,rCT
				double[] subsRates = Dirichlet.sample(instance.mainRand,
						new double[] { 10, 10, 10, 10, 10, 10 });
				// double[] subsRates=new
				// double[]{0.26,0.18,0.17,0.15,0.11,0.13};
				// stationary state frequencies. pi_A, pi_C, pi_G, pi_T
				double[] statFreqs = Dirichlet.sample(instance.mainRand,
						new double[] { 10, 10, 10, 10 });
				// double[] statFreqs=new double[]{0.3,0.2,0.2,0.3};
				double alpha = Sampling.nextDouble(instance.mainRand, 0.1, 0.9); // shape
				// parameter
				// in
				// the
				// Gamma
				// distribution
				if (instance.betterStartVal) {
					double[] subsRatesTmp = new double[subsRates.length];
					for (int i = 0; i < subsRates.length; i++)
						subsRatesTmp[i] = instance.generator.subsRates[i] * 10;
					subsRates = Dirichlet.sample(instance.mainRand,
							subsRatesTmp);
					double[] subStatFreqsTmp = new double[statFreqs.length];
					for (int i = 0; i < statFreqs.length; i++)
						subStatFreqsTmp[i] = instance.generator.stationaryDistribution[i] * 10;
					statFreqs = Dirichlet.sample(instance.mainRand,
							subStatFreqsTmp);
					alpha = Sampling.nextDouble(instance.mainRand,
							Math.max(instance.generator.alpha - 0.1, 0),
							Math.max(instance.generator.alpha + 0.1, 1.0));
				}

				double pInv = 0; // the proportion of invariant sites
				int nCategories = 4;
				int dataRepeatN = nCategories;
				if (pInv > 0)
					dataRepeatN = nCategories + 1;
				MSAPoset align = MSAPoset.parseAlnOrMsfFormats(instance.data);
				Dataset dataset = DatasetUtils.fromAlignment(align,
						instance.sequenceType, dataRepeatN);
				// CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
				// dataset.nSites(), alpha, nCategories, pInv);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				TreeTopologyProcessor trTopo = new TreeTopologyProcessor();
				final int nMCMC = (int) iterScale;
				final int nPMMH = (int) (nMCMC * instance.pmmhPgsExpMix);
				final int nPGS = nMCMC - nPMMH;
				String resultFolder = Execution.getFile("results");
				File output = new File(resultFolder);
				PMMH4GTRIGamma pmmh = new PMMH4GTRIGamma(dataset, options, tdp,
						instance.useTopologyProcessor, trTopo, subsRates,
						statFreqs, alpha, pInv, instance.a_alpha,
						instance.a_pInv, instance.a_statFreqs,
						instance.a_subsRates, nCategories, false,
						instance.isPMCMC4clock, instance.useNJinfo);
				pmmh.setSaveTreesFromPMCMC(instance.saveTreesFromPMCMC);
				if (instance.saveTreesFromPMCMC) {
					String cmdStr = "echo '#NEXUS' > "
							+ instance.nameOfAllTrees;
					String msg0 = IO.call("bash -s", cmdStr, output);
					IO.call("bash -s", "echo 'BEGIN TREES;' >> "
							+ instance.nameOfAllTrees, output);
					IO.call("bash -s", "echo 'TRANSLATE' >> "
							+ instance.nameOfAllTrees, output);
					for (int j = 0; j < dataset.observations().size(); j++)
						IO.call("bash -s", "echo '" + j + " leaf_" + j
								+ ",' >> " + instance.nameOfAllTrees, output);
					IO.call("bash -s",
							"echo ';' >> " + instance.nameOfAllTrees, output);
					pmmh.setNameOfAllTrees(instance.nameOfAllTrees);
				}
				LogInfo.logsForce(" # particles: " + options.nParticles
						+ "; nMCMC:" + nMCMC);
				final int nPMMHburnin = (int) (nPMMH * instance.burninPercent);
				int i = 0;
				while (i < nPMMH) {
					i++;
					System.out.println("PMMH " + i);
					// if (i > nPMMHburnin)
					// pmmh.setProcessTree(true);
					pmmh.next(instance.mainRand);
				}
				RootedTree initTree = pmmh.getRootedTree();
				if (nPMMH == 0)
					initTree = RandomRootedTrees.sampleCoalescent(
							instance.mainRand, align.nTaxa(), 1);
				options.nThreads = instance.nThreads;
				InteractingParticleGibbs4GTRIGamma pg = new InteractingParticleGibbs4GTRIGamma(
						dataset, options, tdp, instance.useTopologyProcessor,
						trTopo, initTree, pmmh.getSubsRates(),
						pmmh.getStateFreqs(), pmmh.getAlpha(), pmmh.getpInv(),
						instance.a_alpha, instance.a_pInv,
						instance.a_statFreqs, instance.a_subsRates,
						nCategories, false, instance.isPMCMC4clock,
						instance.sampleTreeEveryNIter, instance.nCSMC,
						instance.nUCSMC);
				pg.setSaveTreesFromPMCMC(instance.saveTreesFromPMCMC);
				pg.setNameOfAllTrees(instance.nameOfAllTrees);
				i = 0;
				final int nPGSburnin = (int) (nPGS * instance.burninPercent);
				while (i < nPGS) {
					i++;
					System.out.println("Gibbs " + i);
					if (i > nPGSburnin)
						pg.setProcessTree(true);
					pg.next(instance.mainRand);
				}
				if (instance.saveTreesFromPMCMC) {
					IO.call("bash -s", "echo 'END;' >> "
							+ instance.nameOfAllTrees, output);
					String RcmdStr = instance.RcommandDir
							+ "thetaTracePlot.R  \"resultFolder='"
							+ resultFolder + "'\"";
					String msg0 = IO.call("bash -s", RcmdStr, output);
					LogInfo.logs(msg0);
				}

				if (instance.useTopologyProcessor) {
					Counter<UnrootedTree> urtCounter = trTopo.getUrtCounter();
					LogInfo.logsForce("\n Number of unique unrooted trees: "
							+ urtCounter.keySet().size());
					for (UnrootedTree urt : urtCounter.keySet()) {
						LogInfo.logsForce(urt);
						LogInfo.logsForce(urtCounter.getCount(urt));
					}
				}
				return tdp;
			}
		},
		PG {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilter<PartialCoalescentState> pc = new ParticleFilter<PartialCoalescentState>();
				// pc.N = (int) (iterScale * instance.nThousandIters * 1000);
				pc.nThreads = instance.nThreads;
				pc.resampleLastRound = false;
				pc.rand = instance.mainRand;
				pc.verbose = instance.verbose;
				Dataset dataset = DatasetUtils.fromAlignment(instance.data,
						instance.sequenceType);
				// KernelType kernelType = KernelType.PRIOR_PRIOR;
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
				PartialCoalescentState init = PartialCoalescentState
						.initFastState(instance.resampleRoot, dataset, ctmc,
								false);
				NCPriorPriorKernel ppk = new NCPriorPriorKernel(init);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				final double requested = (iterScale * instance.nThousandIters * 1000);
				// pc.N = 1+ (int) Math.pow(requested, instance.pmcmcSMCExpMix);
				pc.N = instance.nParticlesEachStep;
				// pty.smc.PG2 pg = new pty.smc.PG2(pc, ppk, tdp, init);
				pty.smc.PG pg = new pty.smc.PG(pc, ppk, tdp, init);
				final int nMCMC = 1 + (int) Math.pow(requested,
						1.0 - instance.pmcmcSMCExpMix);
				System.out.println("pc.N is " + pc.N + "; nMCMC is " + nMCMC);
				for (int i = 0; i < nMCMC; i++) {
					System.out.println("Iteration " + i);
					pg.next(instance.mainRand);
				}
				// pc.sample(ppk,tdp);
				return tdp;
			}
		},
		/*
		 * BlockPG4GTRIGamma {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilter<PartialCoalescentState> pc = new
		 * ParticleFilter<PartialCoalescentState>(); // pc.N = (int) (iterScale
		 * * instance.nThousandIters * 1000); pc.nThreads = instance.nThreads;
		 * pc.resampleLastRound = false; pc.rand = instance.mainRand; pc.verbose
		 * = instance.verbose;
		 * 
		 * double[] subsRates=instance.generator.subsRates; double[]
		 * statFreqs=instance.generator.stationaryDistribution; double
		 * alpha=instance.generator.alpha; // shape parameter in the Gamma
		 * distribution double pInv=0; // the proportion of invariant sites int
		 * nCategories=4; int dataRepeatN=nCategories;
		 * if(pInv>0)dataRepeatN=nCategories+1; MSAPoset align =
		 * MSAPoset.parseAlnOrMsfFormats(instance.data); Dataset dataset =
		 * DatasetUtils.fromAlignment(align, instance.sequenceType,dataRepeatN);
		 * CTMC ctmc = new CTMC.GTRIGammaCTMC(statFreqs, subsRates, 4,
		 * dataset.nSites(), alpha, nCategories, pInv); PartialCoalescentState
		 * init = PartialCoalescentState.initFastState(dataset, ctmc, false);
		 * NCPriorPriorKernel ppk = new NCPriorPriorKernel(init);
		 * TreeDistancesProcessor tdp = new TreeDistancesProcessor(); final
		 * double requested = (iterScale * instance.nThousandIters * 1000);
		 * //pc.N = 1+ (int) Math.pow(requested, instance.pmcmcSMCExpMix); pc.N
		 * =instance.nParticlesEachStep; //pty.smc.PG pg = new pty.smc.PG(pc,
		 * ppk, tdp, init); pty.smc.PG2 pg = new pty.smc.PG2(pc, ppk, tdp,
		 * init); final int nMCMC = 1 + (int) Math.pow(requested, 1.0 -
		 * instance.pmcmcSMCExpMix);
		 * System.out.println("pc.N is "+pc.N+"; nMCMC is "+nMCMC); for (int i =
		 * 0; i < nMCMC; i++) pg.next(instance.mainRand); // pc.sample(ppk,tdp);
		 * return tdp; } }, NJ {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * NJ nj = new NJ(); MSAPoset align = MSAParser.parseMSA(instance.data);
		 * RootedTree rt = nj.inferTree(align.toMultiAlignmentObject()); final
		 * UnrootedTree urt = UnrootedTree.fromRooted(rt); return new
		 * TreeDistancesProcessor() {
		 * 
		 * @Override public UnrootedTree getConsensus() { return urt; } }; }
		 * 
		 * }, NINJA {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * MSAPoset align = MSAParser.parseMSA(instance.data); final
		 * UnrootedTree tree = Ninja.inferTreeFromNucleotides(align); return new
		 * TreeDistancesProcessor() {
		 * 
		 * @Override public UnrootedTree getConsensus() { return tree; } }; }
		 * 
		 * }, NINJA2 {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * MSAPoset align = MSAParser.parseMSA(instance.data); final
		 * UnrootedTree tree = Ninja.inferTreeFromNucleotides2(align); return
		 * new TreeDistancesProcessor() {
		 * 
		 * @Override public UnrootedTree getConsensus() { return tree; } }; }
		 * 
		 * }, LAZYPMMHNC {
		 * 
		 * @Override public TreeDistancesProcessor doIt(PMCMCExperiments
		 * instance, double iterScale, UnrootedTree goldut, String treeName) {
		 * ParticleFilter<LazyPCS> pc = new ParticleFilter<LazyPCS>(); // pc.N =
		 * (int) (iterScale * instance.nThousandIters * 1000); pc.N=45;
		 * pc.nThreads = 1; // only nThread=1 works properly //
		 * pc.resampleLastRound = false; //
		 * pc.resamplingStrategy=ResamplingStrategy.NEVER; pc.rand =
		 * instance.mainRand;
		 * 
		 * Dataset dataset = DatasetUtils.fromAlignment(instance.data,
		 * instance.sequenceType);
		 * 
		 * // KernelType kernelType = KernelType.PRIOR_PRIOR;
		 * 
		 * double initTrans2tranv= 2.2; TreeDistancesProcessor tdp = new
		 * TreeDistancesProcessor(); LAZYPMMHNC pmmhnc = new LAZYPMMHNC(dataset,
		 * pc, tdp, initTrans2tranv, instance.treeRatePrior);
		 * 
		 * // final double requested = (iterScale * instance.nThousandIters *
		 * 1000); // pc.N = 1+ (int) Math.pow(requested,
		 * instance.pmcmcSMCExpMix); // final int nMCMC = 1 + (int)
		 * Math.pow(requested, 1.0 - instance.pmcmcSMCExpMix); final int
		 * nMCMC=(int) iterScale; System.out.println("# particles: "+pc.N+
		 * "; nMCMC:"+ nMCMC); for (int i = 0; i < nMCMC; i++)
		 * pmmhnc.next(instance.mainRand); return tdp; } };
		 */
		ANNEALING {

			@Override
			public TreeDistancesProcessor doIt(PMCMCExperiments instance,
					double iterScale, UnrootedTree goldut, String treeName) {
				ParticleFilterSMCSampler<UnrootedTreeState> pc = new ParticleFilterSMCSampler<UnrootedTreeState>();
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
				pc.resamplingStrategy = ParticleFilterSMCSampler.ResamplingStrategy.ESS;
				pc.adaptiveTempDiff = instance.adaptiveTempDiff;
				pc.alpha = instance.alphaSMCSampler;

				pc.essRatioThreshold = instance.essRatioThreshold;
				pc.smcSamplerOut = IOUtils.openOutEasy(new File(
						instance.output, "essTempDiff.csv"));

				List<Taxon> leaves = MSAParser.parseMSA(instance.data).taxa();

				RootedTree initTree = TreeGenerators.sampleCoalescent(
						instance.mainRand, leaves, false);
				Gamma exponentialPrior = Gamma.exponential(10.0);
				// StandardNonClockPriorDensity priorDensity = new
				// StandardNonClockPriorDensity(
				// Distrib < Double > branchDistribution);
				StandardNonClockPriorDensity priorDensity = new StandardNonClockPriorDensity(
						exponentialPrior);
				Dataset dataset = DatasetUtils.fromAlignment(instance.data,
						instance.sequenceType);
				// CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());

				System.out.println("instance.csmc_trans2tranv "
						+ instance.csmc_trans2tranv);
				CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites(),
						instance.csmc_trans2tranv);

				UnrootedTreeState ncts = UnrootedTreeState.initFastState(
						UnrootedTree.fromRooted(initTree), dataset, ctmc,
						priorDensity);
				System.out.println(ncts.getLogPrior());

				ProposalDistribution.Options proposalOptions = ProposalDistribution.Util._defaultProposalDistributionOptions;
				LinkedList<ProposalDistribution> proposalDistributions = new LinkedList<ProposalDistribution>();
				ParticleKernel<UnrootedTreeState> ppk = new AnnealingKernel(
						ncts, proposalDistributions, proposalOptions);
				TreeDistancesProcessor tdp = new TreeDistancesProcessor();
				pc.sample(ppk, (ParticleProcessor<UnrootedTreeState>) tdp);
				return tdp;
			}
		};
		abstract TreeDistancesProcessor doIt(PMCMCExperiments instance,
				double iterScale, UnrootedTree goldut, String treeName);
	}

	@SuppressWarnings("unchecked")
	public static <T> double dist(Counter<T> c1, Counter<T> c2) {
		double sum = 0;
		for (T t : union(c1.keySet(), c2.keySet()))
			sum += Math.abs(c1.getCount(t) - c2.getCount(t));
		return sum;
	}

	public static double dist(TreeDistancesProcessor p1,
			TreeDistancesProcessor p2) {
		return dist(p1.getUnrootedCladesPosterior(),
				p2.getUnrootedCladesPosterior());
	}

	@Override
	public void run() {

		if (methods.size() != iterScalings.size())
			throw new RuntimeException(
					"Number of methods and scaling iters should match");

		if (useDataGenerator) {
			LogInfo.logsForce("Generating data...");
			LogInfo.forceSilent = false;
			generator.run();
		} else {

			LogInfo.forceSilent = false;
			File dataDir = new File(Execution.getFile(dataDirName));
			LogInfo.logsForce("Copying data to " + dataDir);
			dataDir.mkdir();
			String str = IO.call("/bin/cp " + dataFile + " " + dataDir + "/"
					+ dataFile.getName());
			LogInfo.logs(str);
		}

		// LogInfo.forceSilent = false;
		output = new File(Execution.getFile("results"));
		output.mkdir();
		treeComparison();

	}

}
