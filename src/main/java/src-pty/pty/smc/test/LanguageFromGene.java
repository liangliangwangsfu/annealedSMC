package pty.smc.test;


import java.util.ArrayList;
import java.util.Random;
import java.util.Set;

import nuts.util.Arbre;
import nuts.util.Counter;
import pty.Train;
import pty.eval.Purity;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.learn.CTMCLoader;
import pty.smc.ConditionalPriorPriorKernel;
import pty.smc.MapLeaves;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.PCSHash;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.BayesRiskMinimizer;
import goblin.Taxon;


public class LanguageFromGene implements Runnable {
//	@Option public String genefile; // use HGDPDataset.path instead
	@Option public String mapfile = "data/language-gene-map.txt";
	@Option public int gibbsIterations = 100;
  @Option public int increaseNSamplesPerGibbsIteration = 0;
	
  private static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
  private static CTMCLoader langParamLoader = new CTMCLoader();
  
  private Dataset langData;
  private MapLeaves ml;
	
	public static double agreementWeight = 1;
	public static final Random rand= new Random(1);
	
	PartialCoalescentState genepcs, langpcs;
	public static void main (String args []){
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
		Execution.run(args, new TestJointModel(),
        "wals", WalsDataset.class, 
        "hgdp", HGDPDataset.class,
        "filter", pf,
        "langparam", langParamLoader,
        "cppk", ConditionalPriorPriorKernel.class);
	}
	

	
	
	public static PartialCoalescentState initCTMCGeneState () { //String filename) {
	  Dataset data = DatasetType.HGDP.loadDataset();
		int nsites = data.nSites();
		double rateScalar = 1e-8;
		double[][] rate = {{-rateScalar,+rateScalar}, {+rateScalar,-rateScalar}} ;
		CTMC ctmc = new CTMC.SimpleCTMC (rate, nsites);
		return  PartialCoalescentState.initState(data,ctmc);
	}
	
	private PartialCoalescentState initLanguageState () {
	  langData = WalsDataset.getPreprocessedCorpus();
	  langParamLoader.setData(langData);
	  CTMC langParam = langParamLoader.load();
	  return PartialCoalescentState.initState(langData, langParam);	
	}
	
	
	public void gibbsSampler (PartialCoalescentState initGeneState, 
			PartialCoalescentState initLangState) {
		// Get an initial Language tree
		ParticleKernel <PartialCoalescentState> ppk = new PriorPriorKernel (initLangState);
    ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> processor =
      SymmetricDiff.createCladeProcessor();
		pf.sample ( ppk , processor);
		Set<Set<Taxon>> currentGeneState = null,
		                   currentLangState = processor.map();
		// initial evaluation
		Arbre<Taxon> reconstruction = Train.outputTree( SymmetricDiff.clades2arbre(processor.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE)), "ling-consensusTree-init");
		if (langData.hasReferenceClusters())
      LogInfo.logs("Purity-init:" + 
          Purity.purity(Arbre.arbre2Tree(reconstruction), 
                        langData.getReferenceClusters()));
		// counters to collect samples across iterations for the MC estimator
		Counter<Set<Set<Taxon>>> allLinguisticSamples = new Counter<Set<Set<Taxon>>>(), 
		                            allGeneSamples = new Counter<Set<Set<Taxon>>>();
 		// do the Gibbs sampling here
		for (int i=0 ; i < gibbsIterations; i++) {
 		  // sample gene
			currentGeneState = sampleBlock(initGeneState, currentLangState, allGeneSamples, false, i);
			// sample languages
			currentLangState = sampleBlock(initLangState, currentGeneState, allLinguisticSamples, true, i);
			pf.N += increaseNSamplesPerGibbsIteration;
		}
	}
	
	/**
	 * Samples one of {bio,ling} node in the graphical model
	 * @param initCoalescentState
	 * @param otherNodeState
	 * @param allSampleForCurrentNode
	 * @param isLang
	 * @param i
	 * @return clade representation of the current node
	 */
	private Set<Set<Taxon>> sampleBlock(PartialCoalescentState initCoalescentState, 
	    Set<Set<Taxon>> otherNodeState, Counter<Set<Set<Taxon>>> allSampleForCurrentNode, boolean isLang, int i)
	{
	  ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> processor 
	    = SymmetricDiff.createCladeProcessor();
	  PCSHash hashProcessor = new PCSHash();
	  ForkedProcessor<PartialCoalescentState> processors = new ForkedProcessor(processor, hashProcessor);
	  ParticleKernel<PartialCoalescentState> pk = new ConditionalPriorPriorKernel(initCoalescentState, otherNodeState, ml, agreementWeight);
    pf.sample( pk, processors);
    allSampleForCurrentNode.incrementAll(processor.getCounter());
    Arbre<Taxon> reconstruction = Train.outputTree(
        SymmetricDiff.clades2arbre(new BayesRiskMinimizer<Set<Set<Taxon>>>(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE).findMin(allSampleForCurrentNode)), 
        (isLang ? "ling" : "bio") + "-consensusTree-" + i);
    if (isLang && langData.hasReferenceClusters())
      LogInfo.logs("Purity-" + i + ":" + 
          Purity.purity(Arbre.arbre2Tree(reconstruction), 
                        langData.getReferenceClusters()));
    LogInfo.logs("Hash-"+ i + (isLang?"-lang":"-bio") + "=" + hashProcessor.getHash());
    return processor.sample(rand);
	}
	
	public void run () {
		PartialCoalescentState geneState = TestBrownianModel.initGeneState (1);
		PartialCoalescentState languageState = initLanguageState ();
		ml = MapLeaves.parse (mapfile);
		gibbsSampler (geneState, languageState);
	}
	
	
}
