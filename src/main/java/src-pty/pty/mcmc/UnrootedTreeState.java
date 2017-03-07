package pty.mcmc;
import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ma.SequenceType;
import monaco.Density;
import nuts.util.CollUtils;
import pty.UnrootedTree;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.WalsDataset.WalsCorpusOperation;
import pty.learn.CTMCLoader;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.models.BrownianModel;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.CTMC;
import pty.smc.models.CTMC.SimpleCTMC;
import pty.smc.models.FastDiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import conifer.Phylogeny;
import conifer.particle.PhyloParticle;
import conifer.particle.PhyloParticleInitContext;
import fig.basic.NumUtils;
import fig.basic.Pair;
import goblin.Taxon;

public final class UnrootedTreeState implements PhyloParticle, Serializable
{
  private static final long serialVersionUID = 1L;

  public UnrootedTreeState()
  {
    this.generationIndex = 0;
  }
  
  
  @Override
  public void init(
      List<Pair<Taxon, LikelihoodModelCalculator>> calculators,
      Density<PhyloParticle> prior,
      PhyloParticleInitContext context)
  {
    if (likelihoodModels != null)
      throw new RuntimeException("Looks like trying to init a UTS constructed from legacy constructors");
    
    likelihoodModels = CollUtils.map();
    for (int i = 0; i < calculators.size(); i++)
      likelihoodModels.put(calculators.get(i).getFirst(), calculators.get(i).getSecond());
    
    // init tree (ORDER IS IMPORTANT--- SHOULD BE CALLED BEFORE PRIOR EVAL)
    this.t = context.getTree().getUnrooted();
    
    this.prior = prior;
    this.logPrior = prior.logDensity(this);

  }

  @Override
  public Phylogeny getPhylogeny()
  {
    return t;
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
  public int generationIndex()
  {
    return generationIndex;
  }
  
  
  
  // What's below can be largely ignored
//  
//  public void changeTree(UnrootedTree t)
//  {
//    this.t = t;
//    if (prior != null)
//      this.logPrior = prior.logDensity(this);
//  }
  
  
  private Density<PhyloParticle> prior;
  private double logPrior;
  private final int generationIndex;
  
  private UnrootedTree t;
//  private final Observations obs;
  private  Map<Taxon,LikelihoodModelCalculator> likelihoodModels;

  private double loglikelihood = Double.NaN;
  
  public UnrootedTree getUnrootedTree() { return t; }
  public UnrootedTree getNonClockTree() { return t; }
  
  public UnrootedTreeState(UnrootedTree t, //Observations obs,
      Map<Taxon, LikelihoodModelCalculator> likelihoodModels,
      Density<PhyloParticle> prior,
      int generationIndex)
  {
    this.generationIndex = generationIndex;
    this.t = t;
//    this.obs = obs;
    this.likelihoodModels = likelihoodModels;
    this.prior = prior;
    if (prior != null)
      this.logPrior = prior.logDensity(this);
  }
  
  @Deprecated
  /**
   * Use init() instead
   * Depracated because prior is not passed along
   */
  public static UnrootedTreeState initFastState(UnrootedTree t, Dataset data, CTMC ctmc)
  {
    List<Taxon> leafNames = new ArrayList<Taxon>();
    Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();
    Map<Taxon, double[][]> observations = data.observations();
    for (Taxon lang : observations.keySet())
    {
      leafNames.add(lang);
      leaves.put(lang, FastDiscreteModelCalculator.observation(ctmc, observations.get(lang), false));
    }
    return new UnrootedTreeState(t, leaves, null, 0);
  }

	@Deprecated
	/**
	 * Use init() instead
	 * Depracated because prior is not passed along
	 *  L. Wang, Nov 19, 2016 
	 */
	public static UnrootedTreeState initFastState(UnrootedTree t, Dataset data,
			CTMC ctmc, Density<PhyloParticle> prior) {
		List<Taxon> leafNames = new ArrayList<Taxon>();
		Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();
		Map<Taxon, double[][]> observations = data.observations();
		for (Taxon lang : observations.keySet()) {
			leafNames.add(lang);
			leaves.put(
					lang,
					FastDiscreteModelCalculator.observation(ctmc,
							observations.get(lang), false));
		}
		return new UnrootedTreeState(t, leaves, prior, 0);
	}

  @Deprecated
  /**
   * Use init() instead
   * Depracated because prior is not passed along
   * added by L. Wang, Sep 24, 2013
   */
  public static UnrootedTreeState initFastState(UnrootedTree t, Map<Taxon, double[][]> observations, CTMC ctmc)
  {
    List<Taxon> leafNames = new ArrayList<Taxon>();
    Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();  
    for (Taxon lang : observations.keySet())
    {
      leafNames.add(lang);
      leaves.put(lang, FastDiscreteModelCalculator.observation(ctmc, observations.get(lang), false));
    }
    return new UnrootedTreeState(t, leaves, null, 0);
  }

  
  public UnrootedTreeState deepClone()
  {
    UnrootedTree nct = new UnrootedTree(t);
    return new UnrootedTreeState(nct, /*this.obs,*/ this.likelihoodModels, this.prior, this.generationIndex);
  }
  
  public UnrootedTreeState copyAndChange(UnrootedTree newTree)
  {
    return new UnrootedTreeState(newTree, /*this.obs,*/ this.likelihoodModels, this.prior, this.generationIndex + 1);
  }
  
	public static UnrootedTreeState fromAlignment(UnrootedTree tree, File align, SequenceType st)
  {
    Dataset data = Dataset.DatasetUtils.fromAlignment(align, st);
    CTMC ctmc = null;
    if (st == SequenceType.RNA || 
        st == SequenceType.DNA)
      ctmc = SimpleCTMC.dnaCTMC(data.nSites());
    else if (st == SequenceType.PROTEIN)
      ctmc = SimpleCTMC.proteinCTMC(data.nSites());
    else
      throw new RuntimeException();

    List<Taxon> leafNames = new ArrayList<Taxon>();
    Map<Taxon, LikelihoodModelCalculator> leaves = CollUtils.map();
    Map<Taxon, double[][]> observations = data.observations();
    for (Taxon lang : observations.keySet())
    {
      leafNames.add(lang);
      leaves.put(lang, FastDiscreteModelCalculator.observation(ctmc, observations.get(lang), false));
    }
		return initFastState(tree, data, ctmc);
//    Map inits = DiscreteModelCalculator.getInit(align, st);
//    return new NonClockTreeState(tree, inits);
  }
  
  /**
   * Sets generation Index to zero
   * @param pcs
   * @return
   */
  public static UnrootedTreeState fromPartialCoalescentState(PartialCoalescentState pcs)
  {
    return fromPartialCoalescentState(pcs, 0);
  }
  public static UnrootedTreeState fromPartialCoalescentState(PartialCoalescentState pcs, int generationIndex)
  {
    return new UnrootedTreeState(
        UnrootedTree.fromRooted(pcs.getFullCoalescentState()), 
        //pcs.getObservations(),
        pcs.getLeafLikelihoodModels(),
        pcs.getPriorDensity(), generationIndex);
  }
  
  @Deprecated
  /**
   * Use init() instead
   * Depracated because prior is not passed along
   */
  public static UnrootedTreeState fromBrownianMotion(UnrootedTree nct, Dataset data, BrownianModel bm)
  {
    // hack:
    PartialCoalescentState pcs = PartialCoalescentState.initState(data, bm, false);
    return new UnrootedTreeState(nct, //pcs.getObservations(), 
        pcs.getLeafLikelihoodModels(), null, 0);
  }
  
  @Deprecated
  /**
   * Use init() instead
   * Depracated because prior is not passed along
   */
  public static UnrootedTreeState fromCTMC(UnrootedTree nct, Dataset data, CTMC ctmc)
  {
    // hack:
    PartialCoalescentState pcs = PartialCoalescentState.initState(data, ctmc);
    return new UnrootedTreeState(nct, //pcs.getObservations(), 
        pcs.getLeafLikelihoodModels(), null, 0);
  }
  

//  
//  public static NonClockTreeState fromMSF(UnrootedTree ut, MSAPoset data, CTMC ctmc)
//  {
////    PartialCoalescentState.
//  }
//  


  public double logLikelihood()
  {
    if (!Double.isNaN(loglikelihood))
      return loglikelihood;
    loglikelihood = computeLogLikelihood(t,likelihoodModels);
    return loglikelihood;
  }
  

  public static double computeLogLikelihood(UnrootedTree t, Map<Taxon,LikelihoodModelCalculator> likelihoodModels)
  {

    // reduce each subtree to a single likelihood model
    final Map<Taxon, LikelihoodModelCalculator> fringe = new LinkedHashMap<Taxon, LikelihoodModelCalculator>();
    
    // this is to make sure that the mergers are done in a deterministic order to make things deterministic
    // (different order will have very similar results modulo different floating point intermeditate steps)
    final ArrayList<Taxon> listOfTaxa = CollUtils.list(likelihoodModels.keySet());
    Collections.sort(listOfTaxa);
    for (Taxon taxon : listOfTaxa)
      fringe.put(taxon, likelihoodModels.get(taxon));
    
    // pick a terminal node, this will be the last one added (need special branch length consideration)
    final Taxon aLeaf = listOfTaxa.get(0);
    
    final Set<Taxon> closedList = new HashSet<Taxon>();
    fringe.remove(aLeaf);
    // combine them, returning their likelihood
    while (fringe.size() > 1)
    {
      // pop any one item
      final Pair<Taxon,Taxon> pair = peek(t,fringe,closedList); //fringe.keySet().iterator().next();
      final Taxon l1 = pair.getFirst(), l2 = pair.getSecond();
      final Taxon parent = parent(t,l1,l2);
      // find relevant b.l.s
      final double 
        BL1 = t.branchLength(l1, parent),
        BL2 = t.branchLength(l2, parent);
      // combine them 
      final LikelihoodModelCalculator 
        LM1 = fringe.get(l1),
        LM2 = fringe.get(l2);
      final LikelihoodModelCalculator combined = 
        LM1.combine(
            LM1, LM2,
            BL1, BL2, false);
      
//      {
//        System.out.println("Combined(" + l1 + "," + l2 + ") -> " + combined.logLikelihood());
//      }
      
      fringe.remove(l1);
      fringe.remove(l2);
      fringe.put(parent, combined);
      closedList.add(l1);
      closedList.add(l2);
    }
    // combine last item
    final Taxon last = fringe.keySet().iterator().next();
    final LikelihoodModelCalculator lastLM = fringe.get(last);
    final LikelihoodModelCalculator aLeafLM = likelihoodModels.get(aLeaf);
    final double halfBL = t.branchLength(last, aLeaf)/2.0;
    final double value = lastLM.combine(
        lastLM, aLeafLM,
        halfBL, halfBL, true).logLikelihood();    
    return value;
  }

  @SuppressWarnings("unchecked")
  private static Taxon parent(UnrootedTree t, Taxon l1, Taxon l2)
  {
    Set<Taxon> inter = CollUtils.inter(t.getTopology().nbrs(l1),t.getTopology().nbrs(l2));
    if (inter.size() != 1) throw new RuntimeException();
    return inter.iterator().next();
  }

  private static Pair<Taxon, Taxon> peek(
      UnrootedTree t,
      Map<Taxon, LikelihoodModelCalculator> fringe,
      Set<Taxon> closedList)
  {
    for (Taxon l1 : fringe.keySet())
    {
      // find node not in closed list
      
      Taxon parent = null;
      {
        List<Taxon> nbrs1 = CollUtils.list(t.getTopology().nbrs(l1));
        Collections.sort(nbrs1);
        
        nbhr:for (Taxon n : nbrs1)
          if (!closedList.contains(n))
          {
            parent = n;
            break nbhr;
          }
      }
      // check that parent has another nbhr in the fringe 
      {
        List<Taxon> nbrs2 = CollUtils.list(t.getTopology().nbrs(parent));
        Collections.sort(nbrs2);
        for (Taxon n : nbrs2)
          if (!n.equals(l1) && fringe.keySet().contains(n))
            return Pair.makePair(l1,n);
      }
    }
    throw new RuntimeException();
  }
  
  @Override
  public String toString()
  {
    String result = t.toString() + "\n";
    result += "LogLikelihood: " + logLikelihood() ;//+ "\t" + "Num sites: " + obs.nSites();
    return result;
  }
  
  private static double harmonicCombine(double v1, double v2, double x1, double x2)
  {
    return (x1/v1 + x2/v2)/(1/v1+1/v2);
  }
  private static double harmonic(double v1, double v2)
  {
    return 1/(1/v1+1/v2);
  }
  
  public static double explicit4LeavesLikelihoodFormula()
  {
    double [] 
      x1 = {0.0,0.9,0.5},  // pop0
      x2 = {0.9,1.0,0.0},  // pop3
      
      x3 = {0.5,0.0,1.0},  // pop1
      x4 = {1.0,0.5,0.9};  // pop2
    
    double 
      v1 = 0.10341776,
      v2 = 0.10207381,
      
      v3 = 0.03380679,
      v4 = 0.02367421,
      
      vcenter = 0.09358231;
    
    double sum = 0.0;
    // big expression
    double bigv = harmonic(v1,v2) + harmonic(v3,v4) + vcenter;
    for (int i = 0; i < x1.length; i++)
      sum += BrownianModelCalculator.logNormalDensity(0, 
          harmonicCombine(v1,v2,x1[i],x2[i]) -
          harmonicCombine(v3,v4,x3[i],x4[i]),
          bigv);
    // small expression b/w x1 & x2
    double smallv = v1 + v2; for (int i = 0; i < x1.length; i++) sum += BrownianModelCalculator.logNormalDensity(0,x1[i]-x2[i], smallv);
    // small expression b/w x3 & x4
           smallv = v3 + v4; for (int i = 0; i < x1.length; i++) sum += BrownianModelCalculator.logNormalDensity(0,x3[i]-x4[i], smallv);
    return sum;
  }
  
  public static double explicit3LeavesLikelihoodFormula(UnrootedTreeState state)
  {
    if (state.likelihoodModels.keySet().size() != 3 ) throw new RuntimeException();
    List<Taxon> leaves = new ArrayList<Taxon>();
    Taxon c = null; // center
    List<double []> xs = new ArrayList<double[]>();
    for (Taxon l : state.getNonClockTree().getTopology().vertexSet())
      if (!state.likelihoodModels.keySet().contains(l)) // is leaf?
        c = l;
      else
      {
        BrownianModelCalculator bm = (BrownianModelCalculator) state.likelihoodModels.get(l);
        leaves.add(l);
        xs.add(bm.message);
      }
    final Taxon 
      l1 = leaves.get(0),
      l2 = leaves.get(1),
      l3 = leaves.get(2);
    final double
      v1 = state.getNonClockTree().branchLength(l1,c),
      v2 = state.getNonClockTree().branchLength(l2,c),
      v3 = state.getNonClockTree().branchLength(l3,c);
    final double []
      x1 = xs.get(0),
      x2 = xs.get(1),
      x3 = xs.get(2);
    final double
      d23 = NumUtils.l2DistSquared(x2,x3),
      d13 = NumUtils.l2DistSquared(x1,x3),
      d12 = NumUtils.l2DistSquared(x1,x2);
    final double 
      sumvs = v1*v2 + v1*v3 + v2*v3,
      sumvds = v1*d23 + v2*d13 + v3*d12;
    final double p = x1.length;
    return -p*Math.log(2*Math.PI) - p*Math.log(sumvs)/2.0 - sumvds/sumvs/2.0;
  }
  
  @SuppressWarnings("unchecked")
  public static void main(String [] args)
  {
    HGDPDataset.path = "data/hgdp/hgdp.ie.phylip";
    BrownianModelCalculator.useVarianceTransform = true;
    Dataset data = new HGDPDataset();
    //
    UnrootedTreeState infileState = UnrootedTreeState.fromBrownianMotion(
        UnrootedTree.fromNewick(new File("data/hgdp/hgdp.ie.tree")), 
        data,
        new BrownianModel(data.nSites(), 1.0));
    
    System.out.println("Code:" + infileState.logLikelihood());
//    System.out.println("Formula:" + explicit3LeavesLikelihoodFormula(infileState));
//    System.out.println("Formula for 4:" + explicit4LeavesLikelihoodFormula());
    if (true) return;
    
    // run a pf
    CTMCLoader langParamLoader = new CTMCLoader();
    WalsDataset.languageFamilyRestriction = "Indo-European";
    WalsDataset.languageListRestriction = new ArrayList(Arrays.asList("fre","eng", "ger", "ita", "bur"));
    WalsDataset.preprocessingSteps = new ArrayList<WalsCorpusOperation>(Arrays.asList(
        WalsCorpusOperation.LIST_RESTRICT,
//        WalsCorpusOperation.FAMILY_RESTRICT,
        WalsCorpusOperation.UNDERDOCUMENTED_LANGS,
        WalsCorpusOperation.UNDERUSED_SITES,
        WalsCorpusOperation.BINARIZE,
        WalsCorpusOperation.REMOVE_DEGENERATE_SITES));
    Dataset langData = WalsDataset.getPreprocessedCorpus();
    langParamLoader.setData(langData);
    CTMC langParam = langParamLoader.load();
    PartialCoalescentState init = PartialCoalescentState.initState(langData, langParam);  
    ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
    ParticleFilter.MAPDecoder<PartialCoalescentState> processor 
      = new ParticleFilter.MAPDecoder<PartialCoalescentState>();
    ParticleKernel<PartialCoalescentState> pk = new PriorPriorKernel(init);
    pf.N = 1;
    pf.sample( pk, processor);
    PartialCoalescentState map = processor.map();
    System.out.println(map.getFullCoalescentState().topology().deepToString());
    // check likelihood of one of final states
    double method1 = map.logLikelihood();
    // convert to nonclocktree
    UnrootedTreeState ncts = UnrootedTreeState.fromPartialCoalescentState(map);
    // compared likelihoods
    System.out.println("Method1:" + method1);
    System.out.println("Method2:" + ncts.logLikelihood());
    // as a sanity for toNewick, convert back and forth and check ll is the same:
    System.out.println(ncts.getNonClockTree().toString());
    String newick = ncts.getNonClockTree().toNewick();
    UnrootedTreeState converted = ncts.copyAndChange(UnrootedTree.fromNewick(newick));
    System.out.println("Sanity:" + converted.logLikelihood());
  }
  public Map<Taxon, LikelihoodModelCalculator> getLikelihoodModels()
  {
    return Collections.unmodifiableMap(likelihoodModels);
  }


}
