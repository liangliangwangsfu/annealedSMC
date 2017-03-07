package pty.smc;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import monaco.Density;
import nuts.io.IO;
import nuts.math.Sampling;
import nuts.maxent.SloppyMath;
import nuts.util.Arbre;
import nuts.util.Arbre.ArbreMap;
import nuts.util.CollUtils;
import pty.Observations;
import pty.RootedTree;
import pty.RootedTree.RootingInfo;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.mcmc.PhyloSampler;
import pty.mcmc.ProposalDistribution.StochasticNearestNeighborInterchangeProposal;
import pty.mcmc.UnrootedTreeState;
import pty.smc.models.BrownianModel;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.FastDiscreteModelCalculator;
import pty.smc.models.ForestModelCalculator;
import pty.smc.models.GTRIGammaFastDiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import pty.smc.models.NoLikelihoodModel;
import pty.smc.test.PhyloHash.Phylo;
import conifer.Phylogeny;
import conifer.data.TaxonIndexedData;
import conifer.particle.FiniteGenerationParticle;
import conifer.particle.PhyloParticle;
import conifer.particle.PhyloParticleInitContext;
import conifer.ssm.SSMModelCalculator;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import goblin.Taxon;

/**
 * A partial coalescent is a forest with each node keeping
 * track of its height  and also caching information for 
 * quickly combining  likelihood model scores
 * 
 * In the rooted case, height is simply time to present, in the unrooted case height is 
 * the min path length to a leaf
 * 
 * Subtrees are constructed with the property that the new, merged subtree should
 * always have a height larger than all the other subtrees. This insures that no two 
 * particle sequences can produce the same state, needed for particle filters.
 * 
 * The initial state is a forest of degenerate trees (as many
 * trees as there are species/languages)
 * @author bouchard
 */
public class  PartialCoalescentState implements PhyloParticle, FiniteGenerationParticle
{
	//  public static int global_id = 0;
	//  public final int id;

	public PartialCoalescentState()
	{
		//    IO.warnOnce("Remove this!");
		//    this.id = global_id++;

		this.topHeight = 0.0;
		this.topNode = null;
		this.oldLogLikelihood = 0.0;
		this.isClock = false;
	}

	public PartialCoalescentState(PartialCoalescentState pcs) {
		this.topHeight = pcs.topHeight;
		this.topNode = pcs.topNode;
		this.oldLogLikelihood = pcs.oldLogLikelihood;
		this.isClock = pcs.isClock;
		this._logLikelihood = pcs._logLikelihood;
		this.data = pcs.data;
		this.logPrior = pcs.logPrior;
		this.logWeightCache = pcs.logWeightCache;
		this.nLeaves = pcs.nLeaves;
		this.obs = pcs.obs;
		this.optimalLogWeight = pcs.optimalLogWeight;
		this.priorFct = pcs.priorFct;
		this.roots = pcs.roots;
	}

	public List<Pair<Taxon, LikelihoodModelCalculator>> calculators;
	public PhyloParticleInitContext context;

	@Override
	public void init(
			List<Pair<Taxon, LikelihoodModelCalculator>> calculators,
			Density<PhyloParticle> prior,
			PhyloParticleInitContext context)
	{
		this.calculators = calculators;
		this.context = context;

		if (this.obs != null)
			throw new RuntimeException("Looks like trying to init PCS created from legacy constructors");
		if (context.getTaxonIndexedData() instanceof Observations)
			this.obs = (Observations) context.getTaxonIndexedData();
		this.data = context.getTaxonIndexedData();

		List<Arbre<CoalescentNode>> roots = new ArrayList<Arbre<CoalescentNode>>();
		for (int i = 0; i < calculators.size(); i++)
		{
			Taxon name = calculators.get(i).getFirst();
			LikelihoodModelCalculator calculator = calculators.get(i).getSecond();
			roots.add(Arbre.arbre(new CoalescentNode(Collections.singleton(name), calculator, 0.0, name, 0.0, 0.0)));
		}
		this.roots = roots;
		this.nLeaves = roots.size();

		this.priorFct = prior;
		this.logPrior = priorFct.logDensity(this);
	}
	@Override
	public Phylogeny getPhylogeny()
	{
		return getFullCoalescentState();
	}
	@Override
	public double getLogLikelihood()
	{
		return logLikelihood();
	}
	@Override
	public double getLogPrior()
	{
		return logPrior;
	}
	@Override
	public int nGenerationsLeft()
	{
		return nIterationsLeft();
	}

	private double logPrior = Double.NaN;
	private Density<PhyloParticle> priorFct;



	// What's below can be largely ignored


	@Option public static boolean disableUnrooted = true;
	@Option public static boolean useExtendedLogLL = false;
	@Option public static boolean useStar = false;
	@Option public static boolean useTemperature = false;  


	public static boolean alwaysComputeTopMessage = false; // this optimization does not work well with PG, so may have to be disabled sometimes


	private final boolean isClock;
	/**
	 * Is this guaranteed to be a clock tree?  Note: even if this return false, the tree might still be clock
	 * @return
	 */
	public boolean isClock() { return isClock; }

	private  TaxonIndexedData data;   // new way
	private  Observations obs; // legacy
	//  private int nTaxa;

	List<Arbre<CoalescentNode>> roots;
	//  private final Set<Set<Taxon>> allClades;
	private final double topHeight;
	private final CoalescentNode topNode;
	private int nLeaves;
	private final double oldLogLikelihood;
	// Useful if we want to associate a weight
	// in addition to the likelihood score e.g. based on the tree topology
	private double logWeightCache;

	private double optimalLogWeight; 

	public transient Set<Set<Taxon>> mapped = null;

	public int nIterationsLeft() { return roots.size() - 1; }
	public boolean isFinalState() { return nIterationsLeft() == 0; }
	public boolean isInitialState() { return topNode == null; }
	public int nRoots() { return roots.size(); }
	public double topHeight() { return topHeight; }

	public int nNonTrivialRoots()
	{
		int i =0;
		for (Arbre<CoalescentNode> subt : roots)
			if (!subt.isLeaf())
				i++;
		return i;
	}

	public int peekNNonTrivialRoots(int i1, int i2)
	{
		int result = nNonTrivialRoots();
		if (roots.get(i1).isLeaf() && roots.get(i2).isLeaf()) 
			return result + 1;
		if (roots.get(i1).isLeaf() || roots.get(i2).isLeaf())
			return result;
		return result - 1;
	}

	public void setWeight (double weight) { logWeightCache = weight;}
	public double getWeight () {return logWeightCache;}
	public int getLeaves () {return nLeaves;}

	public Map<Taxon,LikelihoodModelCalculator> getTopLevelModelCalculators()
	{
		Map<Taxon, LikelihoodModelCalculator> result = CollUtils.map();
		for (Arbre<CoalescentNode> subt : roots)
			result.put(subt.getContents().nodeIdentifier, subt.getContents().likelihoodModelCache);
		return result;
	}

	public Map<Taxon,Integer> rootsTaxa()
	{
		Map<Taxon,Integer> result = CollUtils.map();
		int i = 0;
		for (Arbre<CoalescentNode> root : roots)
			result.put(root.getContents().nodeIdentifier, i++);
		return result;
	}

	/**
	 * What we mean by unlabeled, is that internal nodes are all called null
	 * @return
	 */
	public Arbre<Taxon> getUnlabeledArbre()
	{
		if (!isFinalState()) 
			throw new RuntimeException();
		Arbre<Taxon> result = roots.get(0).preOrderMap(new ArbreMap<CoalescentNode, Taxon>() {
			@Override public Taxon map(Arbre<CoalescentNode> currentDomainNode)
			{ return currentDomainNode.getContents().nodeIdentifier; }
		});
		return result;
	}

	public Phylo<Double> getPhylo()
	{
		if (!isFinalState()) throw new RuntimeException();
		return cnode2phylo(roots.get(0));
	}

	private Phylo<Double> cnode2phylo(Arbre<CoalescentNode> arbre)
	{
		Double h = arbre.getContents().height;
		Set<Phylo<Double>> children = new HashSet<Phylo<Double>>();
		for (Arbre<CoalescentNode> c : arbre.getChildren())
			children.add(cnode2phylo(c));
		return new Phylo<Double>(h, children);
	}
	public List<Arbre<CoalescentNode>> getRoots() {
		return roots;

	}

	public LikelihoodModelCalculator getLikelihoodModelCalculator(int j)
	{
		return roots.get(j).getContents().likelihoodModelCache;
	}
	public boolean isReversible() { return roots.get(0).getContents().likelihoodModelCache.isReversible(); }
	public double getHeight(int j)
	{
		return roots.get(j).getContents().height;
	}

	public RootedTree getFullCoalescentState()
	{
		if (!isFinalState()) 
			throw new RuntimeException();
		final Pair<Arbre<Taxon>,Map<Taxon,Double>> pair
		= getArbreAndBranchLengths();
		return new RootedTree() {
			public Arbre<Taxon> topology()
			{
				return pair.getFirst();
			}
			public Map<Taxon, Double> branchLengths()
			{
				return pair.getSecond();
			}
			@Override
			public int nTaxa()
			{
				return pair.getFirst().nLeaves();
			}
			@Override
			public RootedTree getRooted()
			{
				return this;
			}
			@Override
			public UnrootedTree getUnrooted()
			{
				return UnrootedTree.fromRooted(this);
			}
		};
	}

	public RootedTree getSubtree(int i)
	{
		final Pair<Arbre<Taxon>,Map<Taxon,Double>> pair
		= getArbreAndBranchLengths(i);
		return new RootedTree() {
			public Arbre<Taxon> topology()
			{
				return pair.getFirst();
			}
			public Map<Taxon, Double> branchLengths()
			{
				return pair.getSecond();
			}
			@Override
			public int nTaxa()
			{
				return pair.getFirst().nLeaves();
			}
			@Override
			public RootedTree getRooted()
			{
				return this;
			}
			@Override
			public UnrootedTree getUnrooted()
			{
				return UnrootedTree.fromRooted(this);
			}
		};
	}








	public double nonClockLogWeight()
	{
		if (isClock()) throw new RuntimeException();
		double logSum = Double.NEGATIVE_INFINITY;
		for (Arbre<CoalescentNode>  root : roots)
			if (!root.isLeaf())
			{
				if (root.getChildren().size() != 2)
					throw new RuntimeException();
				final double logCur = + root.getChildren().get(0).getContents().likelihoodModelCache.logLikelihood() 
						+ root.getChildren().get(1).getContents().likelihoodModelCache.logLikelihood()  
						- root.getContents().likelihoodModelCache.logLikelihood();
				logSum = SloppyMath.logAdd(logSum, logCur);
			}
		return -logSum;
	}

	public double logLikelihoodRatio() { 
		//System.out.println("oldLogLikelihood is "+oldLogLikelihood); 
		return logLikelihood() - oldLogLikelihood; }

	private transient double _logLikelihood = Double.NaN;
	public double logLikelihood () 
	{ 
		if (!Double.isNaN(_logLikelihood)) return _logLikelihood;
		if (useStar) return starLogLikelihood();
		else if(useTemperature) return temperatureLogLikelihood();
		else
			return stdLogLikelihood();
	}

	private double stdLogLikelihood()
	{
		if (roots.size() == nLeaves) 
			return 0.0;
		//      _logLikelihood = topNode.likelihoodModelCache.logLikelihood();
		_logLikelihood = 0.0;
		CoalescentNode current;
		for (Arbre<CoalescentNode> root : roots)
			//        if (!root.isLeaf())
		{
			current = root.getContents();
			//        if ((current = root.getContents()) != topNode)
			if (useExtendedLogLL)
				_logLikelihood += current.likelihoodModelCache.extendLogLikelihood(
						topHeight - current.height);  // I think this is ok even in the unrooted case...
			else
				_logLikelihood += current.likelihoodModelCache.logLikelihood();
		}
		return _logLikelihood;

		//    if (isInitialState()) return 0.0;
		//    _logLikelihood = topNode.likelihoodModelCache.logLikelihood();
		//    CoalescentNode current;
		//    for (Arbre<CoalescentNode> root : roots)
		//      if ((current = root.getContents()) != topNode)
		//        _logLikelihood += current.likelihoodModelCache.extendLogLikelihood(
		//            topHeight - current.height);  // I think this is ok even in the unrooted case...
		//    return _logLikelihood;
	}


	public double getTemperature()
	{
		//return 1.0/Math.pow(2.0, (roots.size()-1)*0.5);
		//		double num=1000.0; 
		//		return (num-roots.size()+2)/num;		 	
		double lambda=0.6; 
		double temperature=1-Math.exp(-lambda*(nLeaves-roots.size()+1)); 
		if(Math.abs(temperature-1.0)<1e-2) temperature=1.0;  
		return temperature; 
	}


	private double temperatureLogLikelihood()
	{ 		
		double phi=getTemperature(); 	
		if (roots.size() == nLeaves) 
			return 0.0;
		//      _logLikelihood = topNode.likelihoodModelCache.logLikelihood();
		_logLikelihood = 0.0;
		CoalescentNode current;
		for (Arbre<CoalescentNode> root : roots)
			//        if (!root.isLeaf())
		{
			current = root.getContents();
			//        if ((current = root.getContents()) != topNode)
			if (useExtendedLogLL)
				_logLikelihood += current.likelihoodModelCache.extendLogLikelihood(
						topHeight - current.height);  // I think this is ok even in the unrooted case...
			else
				_logLikelihood += current.likelihoodModelCache.logLikelihood();
		}
		//System.out.println(_logLikelihood+" "+phi); 
		_logLikelihood=_logLikelihood*phi; 		
		return _logLikelihood;
	}


	private double starLogLikelihood()
	{
		double max = Double.NEGATIVE_INFINITY;
		double argmax = -1;
		for (double curDelt = 0.005; curDelt < 20; curDelt *= 2)
		{
			double cur = starLogLikelihood(curDelt);
			//      System.out.println("" + curDelt + "\t" + cur);
			if (cur > max)
			{
				max = cur;
				argmax = curDelt;
			}
		}
		return max;
	}

	private double starLogLikelihood(double plusDelta)
	{
		IO.warnOnce("Using star log likelihood...");
		if (isInitialState()) return 0.0;
		List<Double> deltas = CollUtils.list();
		List<DiscreteModelCalculator> calcs = CollUtils.list();

		for (Arbre<CoalescentNode> root : roots)
		{
			deltas.add(topHeight - root.getContents().height + plusDelta);
			calcs.add((DiscreteModelCalculator) root.getContents().likelihoodModelCache);
		}

		final double result = DiscreteModelCalculator.starCombine(deltas, calcs);


		return result;
	}
	/**
	 * only works when discrete
	 * @return
	 */
	public CTMC getCTMC()
	{
		return ((DiscreteModelCalculator) roots.get(0).getContents().likelihoodModelCache).ctmc;
	}

	public CTMC getFastDiscreteModelCalculatorCTMC()
	{
		return ((FastDiscreteModelCalculator) roots.get(0).getContents().likelihoodModelCache).ctmc;
	}


	public boolean isBrownianMotion()
	{
		return roots.get(0).getContents().likelihoodModelCache instanceof BrownianModelCalculator;
	}
	/**
	 * only works when discrete
	 * @param lang
	 * @return
	 */
	public Map<Taxon,double[]> getObservations(int site)
	{
		if (!isFinalState()) throw new RuntimeException();
		Map<Taxon,double[]> result = new HashMap<Taxon,double[]>();
		for (Arbre<CoalescentNode> node : roots.get(0).nodes())
			if (node.isLeaf())
			{
				DiscreteModelCalculator dmc = (DiscreteModelCalculator) (node.getContents().likelihoodModelCache);
				if (!dmc.isMissing(site))
				{
					double [] cacheCopy = dmc.getCacheCopy(site);
					NumUtils.expNormalize(cacheCopy);
					// NB: nodeId not null b/c it's a leaf
					result.put(node.getContents().nodeIdentifier, cacheCopy);
				}
			}
		return result;
	}

	public Pair<Arbre<Taxon>,Map<Taxon,Double>> getArbreAndBranchLengths(int i) 
	{
		final Map<Taxon,Double> bls = new HashMap<Taxon,Double>();
		Arbre<Taxon> result = roots.get(i).postOrderMap(new ArbreMap<CoalescentNode, Taxon>() {
			int nodeId = 0;
			@Override public Taxon map(Arbre<CoalescentNode> currentDomainNode)
			{ 
				Taxon current = currentDomainNode.getContents().nodeIdentifier; 
				if (current == null) current = new Taxon("internal_" + (nodeId++));
				if (!currentDomainNode.isLeaf())
				{
					List<Taxon> childImage = getChildImage();
					bls.put(childImage.get(0), currentDomainNode.getContents().leftBranchLength);
					bls.put(childImage.get(1), currentDomainNode.getContents().rightBranchLength);
				}
				return current;
			}
		});
		return Pair.makePair(result,bls);
	}

	public Pair<Arbre<Taxon>,Map<Taxon,Double>> getArbreAndBranchLengths() 
	{
		if (!isFinalState()) throw new RuntimeException();
		final Map<Taxon,Double> bls = new HashMap<Taxon,Double>();
		Arbre<Taxon> result = roots.get(0).postOrderMap(new ArbreMap<CoalescentNode, Taxon>() {
			int nodeId = 0;
			@Override public Taxon map(Arbre<CoalescentNode> currentDomainNode)
			{ 
				Taxon current = currentDomainNode.getContents().nodeIdentifier; 
				if (current == null) current = new Taxon("internal_" + (nodeId++));
				if (!currentDomainNode.isLeaf())
				{
					List<Taxon> childImage = getChildImage();
					bls.put(childImage.get(0), currentDomainNode.getContents().leftBranchLength);
					bls.put(childImage.get(1), currentDomainNode.getContents().rightBranchLength);
				}
				return current;
			}
		});
		return Pair.makePair(result,bls);
	}

	@Override
	public String toString()
	{
		StringBuilder result = new StringBuilder();
		int i = 0;
		for (Arbre<CoalescentNode> root : roots)
			result.append("Root " + ((i++)+1) + "/" + roots.size() + "(LL=" + root.getContents().likelihoodModelCache.logLikelihood() + "):\n" + 
					root.preOrderMap(new ArbreMap<CoalescentNode, String>() {
						@Override public String map(Arbre<CoalescentNode> currentDomainNode)
						{
							Taxon nodeId = currentDomainNode.getContents().nodeIdentifier;
							double h = currentDomainNode.getContents().height;
							return "" + (nodeId == null ? "<internal>" : nodeId) + "@h=" + h;
						}
					}).deepToString());
		return result.toString();
	}

	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 */
	public static PartialCoalescentState initialState(
			List<LikelihoodModelCalculator> leaves, 
			List<Taxon> leavesNames,
			Observations obs,
			boolean isClock)
	{
		if (leaves.size() != leavesNames.size()) throw new RuntimeException();
		List<Arbre<CoalescentNode>> roots = new ArrayList<Arbre<CoalescentNode>>();
		for (int i = 0; i < leaves.size(); i++)
			roots.add(Arbre.arbre(new CoalescentNode(Collections.singleton(leavesNames.get(i)), leaves.get(i), 0.0, leavesNames.get(i), 0.0, 0.0)));
		return new PartialCoalescentState(roots, 0.0, leaves.size(), null, /*singletons(rooted,leavesNames),*/ 0.0, obs,isClock, null, null, null, null);
	}


	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 */
	public static PartialCoalescentState initState(Dataset data, BrownianModel bm, boolean resampleRoot)
	{
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		Map<Taxon, double[][]> observations = data.observations();
		for (Taxon lang : observations.keySet())
		{
			leafNames.add(lang);
			double [][] cObs = observations.get(lang);
			double [] converted = new double[cObs.length];
			for (int i = 0; i < converted.length;i++)
				converted[i] = cObs[i][0];
			leaves.add( BrownianModelCalculator.observation (converted, bm, resampleRoot));
			//      leaves.add(DiscreteModelCalculator.observation(ctmc, observations.get(lang)));
		}
		return PartialCoalescentState.initialState(leaves, leafNames, data, (disableUnrooted ? true : false));
	}

	public static PartialCoalescentState initState(List leavesNames)
	{
		return initState(leavesNames, (disableUnrooted ? true : false));
	}

	// for prior sample experiments
	public static PartialCoalescentState initState(List leavesNames, boolean clock)
	{
		//    Observations obs = new Observations() {
		//      public Map<Language, double[][]> observations()
		//      {
		//        return null;
		//      }
		//      public int nCharacter(int site)
		//      {
		//        throw new RuntimeException();
		//      }
		//      public int nSites()
		//      {
		//        throw new RuntimeException();
		//      }
		//    };
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		for (int i = 0; i < leavesNames.size(); i++)
		{
			leafNames.add(new Taxon(leavesNames.get(i).toString()));
			leaves.add(new NoLikelihoodModel());
		}
		return PartialCoalescentState.initialState(leaves, leafNames, null, clock);
	}

	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 */
	public static PartialCoalescentState initState(Dataset data, CTMC ctmc)
	{
		return initState(data,ctmc,null,null);
	}

	public static PartialCoalescentState initForestState(Dataset data, CTMC ctmc, double rootHeight, double langInvRate)
	{
		return initState(data,ctmc,rootHeight, langInvRate);
	}
	private static PartialCoalescentState initState(Dataset data, CTMC ctmc, Double rootHeight, Double langInvRate)
	{
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		Map<Taxon, double[][]> observations = data.observations();
		for (Taxon lang : observations.keySet())
		{
			leafNames.add(lang);
			if (rootHeight == null && langInvRate == null)
				leaves.add(DiscreteModelCalculator.observation(ctmc, observations.get(lang)));
			else if (rootHeight != null && langInvRate != null)
				leaves.add(ForestModelCalculator.observation(ctmc, observations.get(lang), rootHeight, langInvRate));
			else throw new RuntimeException();
		}
		return PartialCoalescentState.initialState(leaves, leafNames, data, true);
	}
	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 */
	public static PartialCoalescentState initFastState(Dataset data, CTMC ctmc)
	{
		return initFastState(false, data,ctmc,true);
	}
	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 */
	public static PartialCoalescentState initFastState(Dataset data, CTMC ctmc, boolean isClock)
	{
		return initFastState(false, data, ctmc, isClock);
	}
	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 */
	public static PartialCoalescentState initFastState(boolean resampleRoot, Dataset data, CTMC ctmc)
	{
		return initFastState(resampleRoot, data, ctmc,true);
	}

	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 */
	public static PartialCoalescentState initFastState(boolean resampleRoot, Dataset data, CTMC ctmc, boolean isClock)
	{
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		Map<Taxon, double[][]> observations = data.observations();
		for (Taxon lang : observations.keySet())
		{
			leafNames.add(lang);
			leaves.add(FastDiscreteModelCalculator.observation(ctmc, observations.get(lang),resampleRoot));
		}
		return PartialCoalescentState.initialState(leaves, leafNames, data, isClock);
	}



	public static PartialCoalescentState initGTRIGammaState(Dataset data, CTMC ctmc, boolean isClock, double alpha, int nCategories)
	{	 
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		Map<Taxon, double[][]> observations = data.observations();
		for (Taxon lang : observations.keySet())
		{
			leafNames.add(lang);
			leaves.add(GTRIGammaFastDiscreteModelCalculator.observation(ctmc, observations.get(lang), alpha, nCategories));
		}
		return PartialCoalescentState.initialState(leaves, leafNames, data, isClock);
	}



	public static PartialCoalescentState initGTRIGammaStateAugMtx(Dataset data, CTMC ctmc, boolean isClock)
	{	 
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		Map<Taxon, double[][]> observations = data.observations();
		DiscreteModelCalculator.allowDiscreteModelCalculator=true;
		for (Taxon lang : observations.keySet())
		{
			leafNames.add(lang);
			leaves.add(DiscreteModelCalculator.observation(ctmc, observations.get(lang)));
		}
		return PartialCoalescentState.initialState(leaves, leafNames, data, isClock);
	}




	public static PartialCoalescentState initSSMState(Map<Taxon,List<String>> data)
	{
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		for (Taxon taxon : data.keySet())
		{
			leafNames.add(taxon);
			leaves.add(new SSMModelCalculator(data.get(taxon)));
		}
		return PartialCoalescentState.initialState(
				leaves, 
				leafNames,
				null,
				true);
	}



	private static Set<Set<Taxon>> singletons(boolean rooted, List<Taxon> leavesNames)
	{
		Set<Set<Taxon>> result = new HashSet<Set<Taxon>>();
		for (Taxon lang : leavesNames)
		{
			result.add(Collections.singleton(lang));
			if (!rooted)
			{
				Set<Taxon> compl = new HashSet<Taxon>(leavesNames);
				compl.remove(lang);
				result.add(compl);
			}
		}
		return result;
	}

	private PartialCoalescentState(
			List<Arbre<CoalescentNode>> roots, 
			double topHeight,
			int nLeaves, CoalescentNode topNode, 
			//      Set<Set<Taxon>> clades,
			double oldLogLikelihood,
			Observations obs,
			boolean isClock,
			TaxonIndexedData data,
			Density<PhyloParticle> priorFct,
			List<Pair<Taxon, LikelihoodModelCalculator>> calculators,
			PhyloParticleInitContext context
			)
	{
		//    IO.warnOnce("Remove this!");
		//    this.id = global_id++;

		this.calculators = calculators;
		this.context = context;
		//  these assertions still hold in unrooted
		assert (topHeight == 0.0 && topNode == null) || 
		topNode.height == topHeight;
		this.roots = roots;
		this.topHeight = topHeight;
		this.nLeaves = nLeaves;
		this.topNode = topNode;
		//    this.allClades = clades;
		this.oldLogLikelihood = oldLogLikelihood;
		this.obs = obs;
		this.isClock = isClock;
		this.data = data;

		// new stuff
		this.priorFct = priorFct;
		if (priorFct != null)
			this.logPrior = priorFct.logDensity(this);
	}

	// this need to be changed
	// NOTE: keep static to avoid memory leaks
	public static class CoalescentNode
	{
		public final LikelihoodModelCalculator likelihoodModelCache;
		public final double height, leftBranchLength, rightBranchLength; // time to present
		public final Taxon nodeIdentifier; 
		public final Set<Taxon> rootedClade;

		private CoalescentNode(
				Set<Taxon> rootedClade,
				LikelihoodModelCalculator likelihoodModelCache, 
				double height, 
				Taxon nodeIdentifier,
				double leftBranchLength, double rightBranchLength)
		{
			this.rootedClade = rootedClade;
			this.likelihoodModelCache = likelihoodModelCache;
			this.height = height;
			this.nodeIdentifier = nodeIdentifier;
			this.leftBranchLength = leftBranchLength;
			this.rightBranchLength = rightBranchLength;
		}
		//    private CoalescentNode(
		//        LikelihoodModelCalculator likelihoodModelCache, 
		//        double height, double leftBranchLength, double rightBranchLength) 
		//    { this(likelihoodModelCache, height, null, leftBranchLength, rightBranchLength); }

		public boolean isLeaf () {
			return height == 0.0;
			//         return nodeIdentifier != null;
		}
	}
	public PartialCoalescentState coalesce(int left, int right, 
			double delta, double leftIncrement, double rightIncrement) 
	{
		return coalesce(left, right, delta, leftIncrement, rightIncrement, null);
	}
	public PartialCoalescentState coalesce(int left, int right, 
			double delta, double leftIncrement, double rightIncrement, Taxon newNode) 
	{
		return (PartialCoalescentState) compute(left, right, delta, leftIncrement, rightIncrement,false, newNode);
	}
	public double peekLogLikelihoodRatio(int left, int right, 
			double delta, double leftIncrement, double rightIncrement) 
	{
		return (Double) compute(left, right, delta, leftIncrement, rightIncrement,true, null);
	}
	private Object compute(int left, int right, 
			double delta, double leftIncrement, double rightIncrement, boolean isPeek, Taxon newNode) 
	{
		{
			//      IO.warnOnce("Clean this hack: PCS");
			if (isClock && (leftIncrement != 0.0 || rightIncrement != 0.0))
				throw new RuntimeException();
			if (!isClock && leftIncrement + rightIncrement == 0.0)
				throw new RuntimeException();
		} 
		//    if (leftIncrement != 0.0 && rightIncrement != 0.0)
		//      throw new RuntimeException();

		if (isFinalState()) throw new RuntimeException();

		//    final boolean isRoot = (roots.size() == 2);
		CoalescentNode  node1 = roots.get(left).getContents(),
				node2 = roots.get(right).getContents();

		final double newTopHeight = delta + topHeight;
		final double branch1 = newTopHeight - node1.height + leftIncrement,
				branch2 = newTopHeight - node2.height + rightIncrement;
		//    System.out.println(newTopHeight+", "+node1.height+", "+node2.height); // all 0's
		//    System.out.println("increment: "+leftIncrement+","+rightIncrement);
		//    System.out.println("two branches: "+branch1+", "+branch2);
		if (isPeek)
		{
			double result = node1.likelihoodModelCache.peekCoalescedLogLikelihood(
					node1.likelihoodModelCache,
					node2.likelihoodModelCache,
					branch1, branch2);
			for (Arbre<CoalescentNode> root : roots)
			{
				CoalescentNode current = root.getContents();
				if (current != node1 && current != node2)
				{
					if (useExtendedLogLL)
						result += current.likelihoodModelCache.extendLogLikelihood(
								newTopHeight - current.height);  
					else
						result += current.likelihoodModelCache.logLikelihood();
				}
			}
			if(useTemperature) result=result*getTemperature();  
            return result - this.logLikelihood(); 
		}
		else
		{
			//      avoid to construct the tables if it's the last
			boolean avoidBuildCache = nIterationsLeft() == 1;
			if (alwaysComputeTopMessage)
				avoidBuildCache = false;
			//      System.out.println(avoidBuildCache);

			// TODO: remove this test code

			//      FastDiscreteModelCalculator.test(node1.likelihoodModelCache, node2.likelihoodModelCache);

			//

			return coalesce(left, right, node1.likelihoodModelCache.combine(node1.likelihoodModelCache,node2.likelihoodModelCache,branch1,branch2,avoidBuildCache), 
					newTopHeight, branch1, branch2, newNode);
		}
	}

	@SuppressWarnings("unchecked")
	private PartialCoalescentState coalesce(int left, int right, LikelihoodModelCalculator contents, 
			double newTopHeight, double branch1, double branch2, Taxon newNode)
	{
		List<Arbre<CoalescentNode>> result = new ArrayList<Arbre<CoalescentNode>>();
		// the trees other than left,right stay the same in the next state 
		for (int i = 0; i < roots.size(); i++)
			if (i != left && i != right)
				result.add(roots.get(i));
		// merge the left and right trees and add it to the forest
		Arbre<CoalescentNode> leftArbre = roots.get(left).copy(),
				rightArbre= roots.get(right).copy();
		Set<Taxon> union = CollUtils.set();
		union.addAll(leftArbre.getContents().rootedClade);
		union.addAll(rightArbre.getContents().rootedClade);
		int ndigits = String.valueOf(nLeaves).length();
		Arbre<CoalescentNode> newArbre = Arbre.arbreWithChildren(
				new CoalescentNode(union,contents, newTopHeight, newNode == null ? new Taxon( /* unique node (only time creates collision!): */ "internal_" +String.format("%0"+ndigits+"d", (nLeaves-nRoots())) + "_" + System.currentTimeMillis() ) : newNode, branch1, branch2),
				leftArbre, rightArbre);
		result.add(newArbre);
		return new PartialCoalescentState(
				result,
				newTopHeight,
				nLeaves, 
				newArbre.getContents(),
				//        nextClades(left,right), 
				this.logLikelihood(),
				this.obs,
				this.isClock,
				this.data,
				this.priorFct,
				this.calculators,
				this.context);
	}
	/**
	 * @deprecated
	 */
	public Set<Taxon> mergedClade(int left, int right)
	{
		throw new RuntimeException();
		//    Set<Taxon> newClade = new HashSet<Taxon>();
		//    newClade.addAll(clade(left));
		//    newClade.addAll(clade(right));
		//    return newClade;
	}
	/**
	 * @deprecated
	 */
	public Set<Set<Taxon>> allClades()
	{
		throw new RuntimeException();
		//    return allClades;
	}
	public Set<Taxon> rootedClade(int i)
	{
		return Collections.unmodifiableSet(roots.get(i).getContents().rootedClade);
	}
	public Set<Set<Taxon>> allRootedClades()
	{
		Set<Set<Taxon>> result = CollUtils.set();
		for (Arbre<CoalescentNode> root : roots)
			allRootedClades(root, result);
		return result;
	}
	private static void allRootedClades(Arbre<CoalescentNode> root,
			Set<Set<Taxon>> result)
	{
		result.add(root.getContents().rootedClade);
		for (Arbre<CoalescentNode> child : root.getChildren())
			allRootedClades(child, result);
	}
	/**
	 * @deprecated
	 */
	public Set<Taxon> clade(int i) 
	{
		throw new RuntimeException();
		//    Set<Taxon> result = new HashSet<Taxon>();
		//    for (CoalescentNode n : roots.get(i).leaveContents()) 
		//      result.add(n.nodeIdentifier);
		//    return result;
	}
	//  private Set<Set<Taxon>> nextClades(int left, int right)
	//  {
	//    Set<Set<Taxon>> result = new HashSet<Set<Taxon>>(allClades);
	//    Set<Taxon> mergedClade = mergedClade(left,right);
	//    if (!isClock)
	//    {
	//      Set<Taxon> compl = CollUtils.unionOfCollection(allClades);
	//      compl.removeAll(mergedClade);
	//      result.add(compl);
	//    }
	//    result.add(mergedClade);
	//    return result;
	//  }

	@Deprecated
	public Observations getObservations()
	{
		return obs;
	}
	@Deprecated
	public int nCharacter(int site)
	{
		return getObservations().nCharacter(site);
	}
	@Deprecated
	public int nSites()
	{
		return getObservations().nSites();
	}
	public Map<Taxon,LikelihoodModelCalculator> getLeafLikelihoodModels()
	{
		Map<Taxon, LikelihoodModelCalculator> result = new HashMap<Taxon,LikelihoodModelCalculator>();
		for (Arbre<CoalescentNode> root : roots)
			for (CoalescentNode leaf : root.leaveContents())
				result.put(leaf.nodeIdentifier, leaf.likelihoodModelCache);
		return result;
	}

	public Density<PhyloParticle> getPriorDensity()
	{
		return priorFct;
	}

	@Override
	public int generationIndex()
	{
		return nLeaves - roots.size();
	}


	public int indexOf(Taxon tax)
	{
		//		System.out.println("looking for "+tax.toString());
		for(int i=0;i<roots.size();i++){
			//			System.out.print(roots.get(i).getContents().nodeIdentifier+"; ");
			if(roots.get(i).getContents().nodeIdentifier.equals(tax)) return i;
		}

		return -1;
	}

	public int randomSelectATreeInForest4NNI(Random rand)
	// Randomly select a tree in a forest (it should have at least 4 leaves).
	{			
		int sampledIndex=-1; 
		List<Integer> candidateList=CollUtils.list();  
		int k=0; 
		for(int i=0;i<roots.size();i++)
			if(roots.get(i).nLeaves()>=4) candidateList.add(k++, i);

		if(candidateList.size()>=1)
		{
			//			System.out.println("Candidate size is: "+candidateList.size());
			int index=Sampling.nextInt(rand, 0, candidateList.size());
			//			System.out.println("index is "+index); 
			sampledIndex =  candidateList.get(index); //List<Integer> sampledIndex = Sampling.sampleWithoutReplacement(rand, this.nRoots(), 1);
			//			System.out.println("sampledIndex is "+sampledIndex); 
		}
		return sampledIndex; 
	}

	public  Pair<Boolean, UnrootedTree> nni4ATreeInForest(Random rand, int sampledIndex)
	{
		RootedTree selectRt=this.getSubtree(sampledIndex);		
		//		Arbre<Taxon> topo=selectRt.topology();
		//		System.out.println(topo.deepToString());

		UnrootedTree ut=selectRt.getUnrooted();
		//		List<Arbre<Taxon>> rootChildren=topo.root().getChildren();
		//		UnorderedPair<Taxon,Taxon> rootEdge=new UnorderedPair<Taxon,Taxon>(rootChildren.get(0).getContents(), rootChildren.get(1).getContents());        	
		PhyloSampler  sampler = new PhyloSampler();
		PhyloSampler.Options phyloSamplerOptions = new PhyloSampler.Options();
		phyloSamplerOptions.rand=rand; 
		sampler.setPhyloSamplerOptions(phyloSamplerOptions); 
		sampler.setOutputText(false);    
		Map<Taxon,double[][]> allObs=this.obs.observations();
		Map<Taxon,double[][]> subObs=new HashMap<Taxon,double[][]>();
		List<Taxon> subleaves=ut.leaves();
		for(int i=0;i<ut.nTaxa();i++)
			subObs.put(subleaves.get(i), allObs.get(subleaves.get(i)));        
		UnrootedTreeState ncts = UnrootedTreeState.initFastState(ut, subObs, this.getFastDiscreteModelCalculatorCTMC()); 
		sampler.init(ncts);
		//boolean accept=sampler.sample(rand);  
		boolean accept=sampler.nnisample(rand, null);
		return Pair.makePair(accept, sampler.getCurrentState().getUnrootedTree()); 
	}

	public  Pair<Boolean, Arbre<CoalescentNode>> nni4ATreeInForest2(Random rand, int sampledIndex) // supposed to be more efficient than nni4ATreeInForest 
	{
		Arbre<CoalescentNode> selectArbre=this.roots.get(sampledIndex);
		RootedTree selectRt=this.getSubtree(sampledIndex);		
		//		Arbre<Taxon> topo=selectRt.topology();
		//		System.out.println(topo.deepToString());		
		UnrootedTree ut=selectRt.getUnrooted();
		double currentLoglike=this.getLikelihoodModelCalculator(sampledIndex).logLikelihood();					
		final double  multiplicativeBranchProposalScaling = 2.0;
		//MultiplicativeBranchProposal proposal=new MultiplicativeBranchProposal(multiplicativeBranchProposalScaling,false); 
		//proposal.selectedEdge = selectedEdge;     
		StochasticNearestNeighborInterchangeProposal proposal=new StochasticNearestNeighborInterchangeProposal(false, multiplicativeBranchProposalScaling); 
		proposal.selectedEdge = null;    
		final Pair<UnrootedTree,Double> result = proposal.propose(ut, rand);	    
		if(result == null) 	return Pair.makePair(false, selectArbre);  // might happen e.g. when trying to do nni with 3 leaves
		final double logProposalRatio = result.getSecond();
		UnrootedTree proposedUt=result.getFirst();    		
		UnorderedPair<Taxon,Taxon> randEdge=proposedUt.randomEdge(rand); 				
		RootingInfo rooting = new RootingInfo(randEdge.getFirst(), randEdge.getSecond(), selectArbre.getContents().nodeIdentifier, rand.nextDouble());
		RootedTree rt = proposedUt.reRoot(rooting);
		Arbre<CoalescentNode> newarbre=rebuildArbreCoales(rt.topology(), rt.branchLengths(), selectArbre.getLeaves()); 
		//		List<Arbre<Taxon>> rootChildren=topo.root().getChildren();
		//		UnorderedPair<Taxon,Taxon> rootEdge=new UnorderedPair<Taxon,Taxon>(rootChildren.get(0).getContents(), rootChildren.get(1).getContents());        	
		double proposeLoglike=newarbre.getContents().likelihoodModelCache.logLikelihood();
		double logPriorRatio=0; // TODO: if any of the branches is changed, logPriorRatio needs to be computed. 
		double logRatio=logPriorRatio+logProposalRatio+proposeLoglike-currentLoglike;  
		final double ratio = Math.min(1,Math.exp(logRatio));
		if (Double.isNaN(ratio)) throw new RuntimeException();
		if (rand.nextDouble() < ratio)
			return Pair.makePair(true, newarbre);
		else 
			return Pair.makePair(false, selectArbre); 
	}
	


	public  Pair<Boolean, UnrootedTree> simplenni4ATreeInForest(Random rand, int sampledIndex)
	{
		RootedTree selectRt=this.getSubtree(sampledIndex);		
		Arbre<Taxon> topo=selectRt.topology();
		//System.out.println(topo.deepToString());
		UnrootedTree ut=selectRt.getUnrooted();
		List<Arbre<Taxon>> rootChildren=topo.root().getChildren();
		UnorderedPair<Taxon,Taxon> rootEdge=new UnorderedPair<Taxon,Taxon>(rootChildren.get(0).getContents(), rootChildren.get(1).getContents());        	
		PhyloSampler  sampler = new PhyloSampler();		
		PhyloSampler.Options phyloSamplerOptions = new PhyloSampler.Options();
		phyloSamplerOptions.rand=rand; 
		sampler.setPhyloSamplerOptions(phyloSamplerOptions); 
		sampler.setOutputText(true);    
		Map<Taxon,double[][]> allObs=this.obs.observations();
		Map<Taxon,double[][]> subObs=new HashMap<Taxon,double[][]>();
		List<Taxon> subleaves=ut.leaves();
		for(int i=0;i<ut.nTaxa();i++)
			subObs.put(subleaves.get(i), allObs.get(subleaves.get(i)));        
		UnrootedTreeState ncts = UnrootedTreeState.initFastState(ut, subObs, this.getFastDiscreteModelCalculatorCTMC()); 
		sampler.init(ncts);
		//boolean accept=sampler.sample(rand);  
		boolean accept=sampler.nnisample(rand, rootEdge);
		return Pair.makePair(accept, sampler.getCurrentState().getUnrootedTree()); 
	}



	public Pair<Pair<Arbre<CoalescentNode>,Arbre<CoalescentNode>>,Pair<Arbre<CoalescentNode>,Arbre<CoalescentNode>>> 
	arbreAfterNNI(Arbre<CoalescentNode> selectArbre, Arbre<Taxon> newtopo)
	{
		List<Arbre<CoalescentNode>> children= selectArbre.getChildren(); 
		Arbre<CoalescentNode> childLeft=children.get(0), childRight=children.get(1);

		//    merge the left and right trees and add it to the forest
		//    Arbre<CoalescentNode> leftArbre = children.get(0).copy(),
		//                          rightArbre= children.get(1).copy();

		List<Arbre<CoalescentNode>> grandChildrenLeft=childLeft.getChildren(), grandChildrenRight=childRight.getChildren();
		//Arbre<CoalescentNode> grandChild0=grandChildrenLeft.get(0).copy(), grandChild1=grandChildrenLeft.get(1).copy(),
		//		grandChild2=grandChildrenRight.get(0).copy(), grandChild3=grandChildrenRight.get(1).copy();

		//Arbre<Taxon> newtopo=rt.topology();
		List<Arbre<Taxon>> newRootChildren=newtopo.getChildren();
		List<Arbre<Taxon>>  grandKids0=newRootChildren.get(0).getChildren(), grandKids1=newRootChildren.get(1).getChildren();
		Taxon grandkid0=grandKids0.get(0).getContents(), grandkid1=grandKids0.get(1).getContents(),
				grandkid2=grandKids1.get(0).getContents(), grandkid3=grandKids1.get(1).getContents();

		List<Arbre<CoalescentNode>> arbreList=CollUtils.list(); 
		arbreList.add(grandChildrenLeft.get(0).copy());
		arbreList.add(grandChildrenLeft.get(1).copy());
		arbreList.add(grandChildrenRight.get(0).copy());
		arbreList.add(grandChildrenRight.get(1).copy()); 

		int[] indice0=new int[2], indice1=new int[2]; 
		for(int i=0;i<arbreList.size();i++)
		{
			if(grandkid0==arbreList.get(i).getContents().nodeIdentifier) indice0[0]=i;
			else
				if(grandkid1==arbreList.get(i).getContents().nodeIdentifier) indice0[1]=i;
				else
					if(grandkid2==arbreList.get(i).getContents().nodeIdentifier) indice1[0]=i;
					else indice1[1]=i;     	
		}		
		return Pair.makePair(Pair.makePair(arbreList.get(indice0[0]),arbreList.get(indice0[1])), Pair.makePair(arbreList.get(indice1[0]),arbreList.get(indice1[1])));
	}


	public int randomSelectATreeInForest4ASimpleNNI(Random rand)
	// Randomly select a tree in a forest (it should have at least 4 leaves).
	{			
		int sampledIndex=-1; 
		List<Integer> candidateList=CollUtils.list();  
		int k=0; 
		for(int i=0;i<roots.size();i++)
		{
			if(!roots.get(i).isLeaf()){
				List<Arbre<CoalescentNode>> children=roots.get(i).getChildren();
				if(!children.get(0).isLeaf() && !children.get(1).isLeaf()) 	candidateList.add(k++, i);
			}
		}		
		if(candidateList.size()>=1)
		{
			System.out.println("Candidate size is: "+candidateList.size());
			int index=Sampling.nextInt(rand, 0, candidateList.size());
			System.out.println("index is "+index); 
			sampledIndex =  candidateList.get(index); //List<Integer> sampledIndex = Sampling.sampleWithoutReplacement(rand, this.nRoots(), 1);
			System.out.println("sampledIndex is "+sampledIndex); 
		}
		return sampledIndex; 
	}

	public  Arbre<CoalescentNode> rebuildArbreCoales(Arbre<Taxon> topo, Map<Taxon,Double> rtbr, List<Arbre<CoalescentNode>> leaves)  //TODO: make it more efficient
	{
		Map<Taxon, Arbre<CoalescentNode>> leavesMap=CollUtils.map();		
		for(int i=0;i<leaves.size();i++)		
			leavesMap.put(leaves.get(i).getContents().nodeIdentifier, leaves.get(i));  		

		List<Arbre<Taxon>> children=topo.getChildren();
		Arbre<Taxon> node1=children.get(0), node2=children.get(1); 
		Taxon tax1=node1.getContents(), tax2=node2.getContents();
		Arbre<CoalescentNode> arbre1, arbre2; 
		if(node1.isLeaf()) arbre1=leavesMap.get(tax1);
		else arbre1=rebuildArbreCoales(node1,  rtbr, leaves);

		if(node2.isLeaf()) arbre2=leavesMap.get(tax2);
		else arbre2=rebuildArbreCoales(node2,  rtbr, leaves);

		return comineArbres(arbre1, arbre2, topo.getContents(), rtbr.get(tax1), rtbr.get(tax2), false);
	}

	public PartialCoalescentState oneMcmcStep(Random rand)
	{				
		int sampledIndex=randomSelectATreeInForest4NNI(rand);
		if(sampledIndex==-1) return this;
		Arbre<CoalescentNode> selectArbre=this.roots.get(sampledIndex);					 						       

		//		List<Arbre<CoalescentNode>> arbreInPostOrder=selectArbre.nodesInPostOrder();
		Pair<Boolean, Arbre<CoalescentNode>> nniRe=nni4ATreeInForest2(rand, sampledIndex); 		

		//		RootedTree selectRt=this.getSubtree(sampledIndex);		
		//		Arbre<Taxon> topo=selectRt.topology();
		//		System.out.println(topo.deepToString());
		//		System.out.println(proposedUt);

		List<Arbre<CoalescentNode>> newroots = new ArrayList<Arbre<CoalescentNode>>();
		for(int i=0;i<roots.size();i++)
		{
			if(i!=sampledIndex) newroots.add(this.roots.get(i)); 
			else newroots.add(nniRe.getSecond());
		}

		return new PartialCoalescentState(
				newroots,
				0,
				nLeaves, 
				selectArbre.getContents(),
				//	        nextClades(left,right), 
				this.logLikelihood(),
				this.obs,
				this.isClock,
				this.data,
				this.priorFct,
				this.calculators,
				this.context);	
	}



	public PartialCoalescentState mcmcStep(Random rand)
	{				
		int sampledIndex=randomSelectATreeInForest4ASimpleNNI(rand);
		if(sampledIndex==-1) return this;
		Arbre<CoalescentNode> selectArbre=this.roots.get(sampledIndex);

		//List<Arbre<CoalescentNode>> arbreInPostOrder=selectArbre.nodesInPostOrder();


		//System.out.println(this.getSubtree(sampledIndex).topology().deepToString()); 

		Pair<Boolean, UnrootedTree> nniRe=simplenni4ATreeInForest(rand, sampledIndex); 

		if(!nniRe.getFirst()) return this; 

		UnrootedTree proposedUt=nniRe.getSecond();    
		//System.out.println(proposedUt);				

		List<Arbre<CoalescentNode>> rootChildren=selectArbre.getChildren();
		CoalescentNode node1=rootChildren.get(0).getContents(), node2=rootChildren.get(1).getContents(); 
		Taxon tax1=node1.nodeIdentifier, tax2=node2.nodeIdentifier;

		RootingInfo rooting = new RootingInfo(tax1, tax2, selectArbre.getContents().nodeIdentifier, rand.nextDouble());
		RootedTree rt = proposedUt.reRoot(rooting);
		//System.out.println(rt);

		Pair<Pair<Arbre<CoalescentNode>,Arbre<CoalescentNode>>,Pair<Arbre<CoalescentNode>,Arbre<CoalescentNode>>> arbres=arbreAfterNNI(selectArbre, rt.topology());

		Map<Taxon,Double> rtbr=rt.branchLengths();

		Arbre<CoalescentNode> arbre11=arbres.getFirst().getFirst(),arbre12=arbres.getFirst().getSecond(),
				arbre21=arbres.getSecond().getFirst(), arbre22=arbres.getSecond().getSecond(); 

		double  branch11 = rtbr.get(arbre11.getContents().nodeIdentifier), branch12 = rtbr.get(arbre12.getContents().nodeIdentifier),
				branch21 = rtbr.get(arbre21.getContents().nodeIdentifier), branch22 = rtbr.get(arbre22.getContents().nodeIdentifier);
		Arbre<CoalescentNode> arbre1=comineArbres(arbre11, arbre12, tax1, branch11, branch12, false), 
				arbre2=comineArbres(arbre21, arbre22, tax2, branch21, branch22, false);

		Arbre<CoalescentNode> newarbre=comineArbres(arbre1, arbre2, selectArbre.getContents().nodeIdentifier, rtbr.get(tax1), rtbr.get(tax2), false);

		List<Arbre<CoalescentNode>> newroots = new ArrayList<Arbre<CoalescentNode>>();
		for(int i=0;i<roots.size();i++)
		{
			if(i!=sampledIndex) newroots.add(this.roots.get(i)); 
			else newroots.add(newarbre);
		}

		return new PartialCoalescentState(
				newroots,
				0,
				nLeaves, 
				selectArbre.getContents(),
				//	        nextClades(left,right), 
				this.logLikelihood(),
				this.obs,
				this.isClock,
				this.data,
				this.priorFct,
				this.calculators,
				this.context);	
	}


	public Arbre<CoalescentNode> comineArbres(Arbre<CoalescentNode> arbre1, Arbre<CoalescentNode> arbre2, Taxon newTopTax, double branch1, double branch2, boolean copyArbres)
	{		
		Arbre<CoalescentNode> leftArbre=arbre1, rightArbre=arbre2; 
		if(copyArbres)
		{
			leftArbre=arbre1.copy(); 
			rightArbre=arbre2.copy();
		}	

		Set<Taxon> union = CollUtils.set();
		union.addAll(leftArbre.getContents().rootedClade);
		union.addAll(rightArbre.getContents().rootedClade);
		//int ndigits = String.valueOf(nLeaves).length();    

		LikelihoodModelCalculator likeLeft=leftArbre.getContents().likelihoodModelCache, likeRight=rightArbre.getContents().likelihoodModelCache; 

		boolean avoidBuildCache=false;   //TODO: check if this is false or true!!

		LikelihoodModelCalculator contents=likeLeft.combine(likeLeft,likeRight,branch1,branch2,avoidBuildCache);
		double newTopHeight=0; 
		return Arbre.arbreWithChildren(
				new CoalescentNode(union,contents, newTopHeight,  newTopTax, branch1, branch2),
				leftArbre, rightArbre);
	}

}
