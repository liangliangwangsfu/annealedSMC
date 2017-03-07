package pty.smc.test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import ma.newick.NewickParser;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.Tree;
import nuts.util.Arbre.ArbreMap;
import pty.RootedTree;
import pty.Train;
import pty.eval.Purity;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.LeaveOneOut;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.learn.CTMCLoader;
import pty.smc.ConditionalPriorPriorKernel;
import pty.smc.MapLeaves;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.PCSHash;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.models.BrownianModel;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import pty.smc.models.ProductModel;
import pty.smc.test.TestBrownianModel.KernelType;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.BayesRiskMinimizer;
import goblin.Taxon;


public class TestJointModel implements Runnable {
//	@Option public String genefile; // use HGDPDataset.path instead
	@Option public String mapfile = "data/language-gene-map.txt";
	@Option public int gibbsIterations = 100;
  @Option public int increaseNSamplesPerGibbsIteration = 0;
  @Option public double variance = 0.1;
  @Option public boolean testAgainstFixedTree = false;
  @Option public String fixedTreePath = "data/hgdp/contml.all.newick";
  @Option public double agreementWeight = 1;
  @Option public boolean softAgreement = true;
	
  private static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
  private static CTMCLoader langParamLoader = new CTMCLoader();
  
  private Dataset langData;
  private MapLeaves ml;
	private CTMC langParam;
	
	public static final Random rand= new Random(1);
	
	PartialCoalescentState genepcs, langpcs;
	public static void main (String args []){
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
		Execution.run(args, new TestJointModel(),
        "wals", WalsDataset.class, 
        "hddp", HGDPDataset.class,
        "filter", pf,
        "langparam", langParamLoader,
        "cppk", ConditionalPriorPriorKernel.class);
	}
	


	
	private PartialCoalescentState initLanguageState () {
	  langData = WalsDataset.getPreprocessedCorpus();
	  langParamLoader.setData(langData);
	  langParam = langParamLoader.load();
	  return PartialCoalescentState.initState(langData, langParam);	
	}
	
	
	public void gibbsSampler (PartialCoalescentState initGeneState, 
			PartialCoalescentState initLangState) {
		// Get an initial Language tree
		ParticleKernel <PartialCoalescentState> ppk = (ConditionalPriorPriorKernel.usesPriorPost ? 
		    KernelType.PRIOR_POST2.load(initLangState,null) : 
		    KernelType.PRIOR_PRIOR.load(initLangState,null));
		    //new PriorPriorKernel (initLangState);
    ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> processor =
      SymmetricDiff.createCladeProcessor();
		pf.sample ( ppk , processor);
		Set<Set<Taxon>> currentGeneState = null,
		                   currentLangState = processor.map();
		// initial evaluation
		Arbre<Taxon> reconstruction = Train.outputTree( SymmetricDiff.clades2arbre(processor.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE)), "ling-consensusTree-init");
		if (langData.hasReferenceClusters())
		{
		  Tree<Taxon> recon = Arbre.arbre2Tree(reconstruction);
		  Map<Taxon,String> allLabels = langData.getReferenceClusters();
		  LogInfo.logs("Labels used for evaluation:" + Purity.partitionsUsedForEval(recon,allLabels)); 
      LogInfo.logs("Purity-init:" +  Purity.purity(recon, allLabels));
		}
		// counters to collect samples across iterations for the MC estimator
		Counter<Set<Set<Taxon>>> allLinguisticSamples = new Counter<Set<Set<Taxon>>>(), 
		                            allGeneSamples = new Counter<Set<Set<Taxon>>>();
 		// do the Gibbs sampling here
		for (int i=0 ; i < gibbsIterations; i++) {
 		  // sample gene
		  if (testAgainstFixedTree)
		    currentGeneState = loadFixedGeneState(fixedTreePath);
		  else
		    currentGeneState = sampleBlock(initGeneState, currentLangState, allGeneSamples, false, i);
			// sample languages
			currentLangState = sampleBlock(initLangState, currentGeneState, allLinguisticSamples, true, i);
			pf.N += increaseNSamplesPerGibbsIteration;
		}
	}
	
	private Set<Set<Taxon>> _fixedGS = null;
	private Set<Set<Taxon>> loadFixedGeneState(String fixedTreePath)
  {
    if (_fixedGS != null) return _fixedGS;
    try
    {
      NewickParser np = new NewickParser(IOUtils.openIn(fixedTreePath));
      Tree<String> tree = np.parse();
      _fixedGS = SymmetricDiff.cladesFromUnrooted(
          Arbre.tree2Arbre(tree).postOrderMap(
              new ArbreMap<String,Taxon>() {
                public Taxon map(Arbre<String> d) { return new Taxon(d.getContents()); }}));
      LogInfo.logs("Fixed constraints:" + _fixedGS);
      return _fixedGS;
    } catch (Exception e) { throw new RuntimeException(e); }
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
	  ParticleFilter.MAPDecoder<PartialCoalescentState> mapDecoder = new ParticleFilter.MAPDecoder<PartialCoalescentState>();
	  ForkedProcessor<PartialCoalescentState> processors = new ForkedProcessor(processor, hashProcessor, mapDecoder);
	  ParticleKernel<PartialCoalescentState> pk = new ConditionalPriorPriorKernel(initCoalescentState, otherNodeState, ml, agreementWeight);
    pf.sample( pk, processors);
    allSampleForCurrentNode.incrementAll(processor.getCounter());
    RootedTree map = mapDecoder.map().getFullCoalescentState();
    final String prefix = (isLang ? "ling" : "bio");
    Train.outputTree (map.topology(), prefix + "-mapTree-" + i, map.branchLengths());
    Arbre<Taxon> reconstruction = Train.outputTree(
        SymmetricDiff.clades2arbre(new BayesRiskMinimizer<Set<Set<Taxon>>>(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE).findMin(allSampleForCurrentNode)), 
        prefix + "-consensusTree-" + i);
    if (isLang && langData.hasReferenceClusters())
      LogInfo.logs("Purity-" + i + ":" + 
          Purity.purity(Arbre.arbre2Tree(reconstruction), 
                        langData.getReferenceClusters()));
    if (isLang)
      LogInfo.logs("LOO-" + i + ":" + LeaveOneOut.loo(mapDecoder.map()));
    LogInfo.logs("Hash-"+ i + (isLang?"-lang":"-bio") + "=" + hashProcessor.getHash());
    LogInfo.logs("TopologyDistributionLargestPrt=" + processor.getCounter().max());
    return processor.sample(rand);
	}
	
	public void run () {
		PartialCoalescentState geneState = TestBrownianModel.initGeneState (variance);
		PartialCoalescentState languageState = initLanguageState ();
		ml = MapLeaves.parse (mapfile);
		if (softAgreement)
		  gibbsSampler (geneState, languageState);
		else
		  hardAgreement();
	}

  private void hardAgreement()
  {
    PartialCoalescentState jointInit = initHardState();
    ParticleKernel <PartialCoalescentState> ppk = (ConditionalPriorPriorKernel.usesPriorPost ? 
        KernelType.PRIOR_POST2.load(jointInit,null) : 
        KernelType.PRIOR_PRIOR.load(jointInit,null));
    // processors
    ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> processor 
      = SymmetricDiff.createCladeProcessor();
    PCSHash hashProcessor = new PCSHash();
    ParticleFilter.MAPDecoder<PartialCoalescentState> mapDecoder = new ParticleFilter.MAPDecoder<PartialCoalescentState>();
    ForkedProcessor<PartialCoalescentState> processors = new ForkedProcessor(processor, hashProcessor, mapDecoder);
    // sample!
    pf.sample( ppk, processors);
    // log
    RootedTree map = mapDecoder.map().getFullCoalescentState();
    Train.outputTree (map.topology(), "mapTree", map.branchLengths());
    Arbre<Taxon> reconstruction = Train.outputTree(
        SymmetricDiff.clades2arbre(processor.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE)),
        "consensusTree");
    if (langData.hasReferenceClusters())
      LogInfo.logs("Purity:" + 
          Purity.purity(Arbre.arbre2Tree(reconstruction), 
                        langData.getReferenceClusters()));
    LogInfo.logs("Hash=" + hashProcessor.getHash());
  }
	
  @Deprecated
  /**
   * Use init() instead
   * Depracated because prior is not passed along
   */
  private PartialCoalescentState initHardState()
  {
    List<Taxon> leafNames = new ArrayList<Taxon>();
    List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
    Map<Taxon, double[][]> 
      langObservations = langData.observations(),
      bioObservations  = DatasetType.HGDP.loadDataset().observations();
    if (langObservations.keySet().size() != bioObservations.keySet().size())
      throw new RuntimeException();
    // will need the correspondance file...
    for (final Taxon lang : langObservations.keySet())
    {
      final Taxon bioEq = ml.translate(lang);
      leafNames.add(lang);
      List<LikelihoodModelCalculator> models = CollUtils.list();
      // create the discrete (ling) leaf
      models.add(DiscreteModelCalculator.observation(langParam, langObservations.get(lang)));
      // create the freq (bio) leaf
      double [][] cObs = bioObservations.get(bioEq);
      double [] converted = new double[cObs.length];
      for (int i = 0; i < converted.length;i++)
        converted[i] = cObs[i][0];
      BrownianModel bm = new BrownianModel(converted.length, variance);
      models.add(BrownianModelCalculator.observation(converted, bm, false));
      leaves.add(new ProductModel(models));
    }
    return PartialCoalescentState.initialState(leaves, leafNames, null, true);
  }
	
}
