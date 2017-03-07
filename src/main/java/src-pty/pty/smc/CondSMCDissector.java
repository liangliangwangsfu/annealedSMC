package pty.smc;
import fig.basic.NumUtils;
import fig.basic.Pair;
import fig.basic.UnorderedPair;
import fig.prob.SampleUtils;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import pty.Observations;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.RootedTree.RootingInfo;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.ParticleFilter.StoreProcessor;
import pty.smc.models.LikelihoodModelCalculator;
import nuts.math.Sampling;
import nuts.util.CollUtils.*;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.Tree;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class CondSMCDissector
{
  private final UnrootedTree originalTree;
  private final Taxon leaf1, leaf2;
  private List<Arbre<Taxon>> backbone;
  private int anchorIndex = -1;
  private Taxon anchor;
  private Set<Taxon> sampled = null, conditioned = null;
  private Arbre<Taxon> rootedAtAnchor;
  
  public CondSMCDissector(UnrootedTree originalTree, Taxon leaf1, Taxon leaf2)
  {
    this.originalTree = originalTree;
    this.leaf1 = leaf1;
    this.leaf2 = leaf2;
    findBackbone();
  }

  private void findBackbone()
  {
    if (leaf1 == leaf2)
      throw new RuntimeException();
    Tree<Taxon> t = originalTree.toTree(leaf1);
    Arbre<Taxon> iter = Arbre.findFirstNodeWithContents(Arbre.tree2Arbre(t), leaf2);
    backbone = list();
    backbone.add(iter);
    while (!iter.isRoot())
    {
      iter = iter.getParent();
      backbone.add(iter);
    }
  }
  
  private Arbre<Taxon> anchorNeighbor(int i)
  {
    if (i == 0 || i == backbone.size() - 1)
      throw new RuntimeException();
    for (Arbre<Taxon> child : backbone.get(i).getChildren())
      if (child != backbone.get(i-1))
        return child;
    throw new RuntimeException();
  }
  
  public boolean sampleAnchor(Random rand)
  {
    if (anchorIndex != -1)
      throw new RuntimeException();
    final int len = backbone.size();
    double [] prs = new double[len]; 
    for (int i = 1; i < len - 1; i++)
      prs[i] = anchorNeighbor(i).nodes().size() - 1;
    boolean success = NumUtils.normalize(prs);
    if (!success)
      return false;
    anchorIndex = SampleUtils.sampleMultinomial(rand, prs);
    anchor = backbone.get(anchorIndex).getContents();
    //
    sampled = set(anchorNeighbor(anchorIndex).nodeContents());
    conditioned =  set(backbone.get(backbone.size()-1).nodeContents());
    conditioned.removeAll(sampled);
    sampled.add(anchor); // becomes a pseudo-leaf
    sampled.add(Taxon.dummy);
    rootedAtAnchor = Arbre.tree2Arbre(originalTree.toTree(anchor));
    return true;
  }
  
  private Arbre<Taxon> sampleRooting(Random rand)
  {
    List<UnorderedPair<Taxon,Taxon>> candidates = list();
    for (UnorderedPair<Taxon,Taxon> edge : originalTree.edges())
      if (sampled.contains(edge.getFirst()) && sampled.contains(edge.getSecond()))
        candidates.add(edge);
    int index = rand.nextInt(candidates.size());
    RootingInfo rooting = new RootingInfo(candidates.get(index).getFirst(), candidates.get(index).getSecond(), Taxon.dummy, 0.5);
    RootedTree rt = originalTree.reRoot(rooting);
    return rt.topology();
  }
  
  private LikelihoodModelCalculator fixedModelCalculator(PartialCoalescentState model)
  {
    GuidedCoalescentCreator gcc = new GuidedCoalescentCreator(conditioned, rootedAtAnchor, model, new Random(1) /* does not actually depend on it*/,
      originalTree, null, null);
    List<PartialCoalescentState> list = gcc.compute();
    return list.get(list.size() -1 ).getLikelihoodModelCalculator(0);
  }
  
  public List<PartialCoalescentState> currentSampledPath(PartialCoalescentState model, Random rand)
  {
    // first, compute message from conditioning set going out of anchor
    LikelihoodModelCalculator addit = fixedModelCalculator(model);
    return new GuidedCoalescentCreator(
        sampled,
        sampleRooting(rand),
        model, 
        rand,
        originalTree,
        anchor,
        addit).compute();
    
  }
  
  public UnrootedTree reconstitute(PartialCoalescentState _sampled)
  {
    UnrootedTree sampled = UnrootedTree.fromRooted(_sampled.getFullCoalescentState());
    
    // NB: depends on having unique internal node names (PartialCoalescent modified accordingly)
    Set<Taxon> inSampled = set();
    for (Taxon t : sampled.getTopology().vertexSet())
      if (sampled.getTopology().nbrs(t).size() != 1)
        inSampled.add(t);
    if (CollUtils.intersects(inSampled, conditioned))
      throw new RuntimeException();
    
    // start with topology
    Arbre<Taxon> reRooted = Arbre.tree2Arbre(sampled.toTree(anchor));
    if (reRooted.getChildren().size() != 1)
      throw new RuntimeException();
    Arbre<Taxon> stemmedRerooted = reRooted.getChildren().get(0);
    List<Arbre<Taxon>> newChildren = list();
    if (rootedAtAnchor.getChildren().size() != 3)
      throw new RuntimeException();
    for (Arbre<Taxon> oldChildren : rootedAtAnchor.getChildren())
      if (this.sampled.contains(oldChildren.getContents()))
        newChildren.add(stemmedRerooted.copy());
      else if (this.conditioned.contains(oldChildren.getContents()))
        newChildren.add(oldChildren.copy());
      else
        throw new RuntimeException();
    Arbre<Taxon> newTree = Arbre.arbre(anchor, newChildren);
    
    // branch lengths now...
    Map<Taxon,Double> bls = map();
    for (Arbre<Taxon> subt : rootedAtAnchor.nodes())
      if (!subt.isRoot())
        if ( this.conditioned.contains(   subt.getContents()))
          bls.put(subt.getContents(), originalTree.branchLength(subt.getContents(), subt.getParent().getContents()));
    for (Arbre<Taxon> subt : reRooted.nodes())
      if (!subt.isRoot())
        bls.put(subt.getContents(), sampled.branchLength(subt.getContents(), subt.getParent().getContents()));
    
    RootedTree resultRT = RootedTree.Util.create(newTree, bls);
    return UnrootedTree.fromRooted(resultRT);
  }
  
  private static class GuidedCoalescentCreator
  {
    private final List<PartialCoalescentState> list = list();
    private final Set<Taxon> taxaToConsider;
    private final Arbre<Taxon> topology;
    private final PartialCoalescentState model;
    private final Random rand;
    private final UnrootedTree originalTree;
    
    private final Taxon additionalLeaf;
    private final LikelihoodModelCalculator additionalLLC;
    
    private PartialCoalescentState current = null;
    private double halfTopBranch = -1;
    
    private GuidedCoalescentCreator(
        Set<Taxon> taxaToConsider,
        Arbre<Taxon> topology, 
        PartialCoalescentState model, 
        Random rand,
        UnrootedTree originalTree, 
        Taxon additionalLeaf,
        LikelihoodModelCalculator additionalLLC)
    {
      this.taxaToConsider = taxaToConsider;
      this.topology = topology;
      this.model = model;
      this.rand = rand;
      this.originalTree = originalTree;
      this.additionalLeaf = additionalLeaf;
      this.additionalLLC = additionalLLC;
    }

    private List<PartialCoalescentState> compute()
    {
      init();
      while (!current.isFinalState())
        iterate();
      return list;
    }
    
    private void init()
    {
      if (model.isClock())
        throw new RuntimeException();
      if (topology.getChildren().size() == 2)
        halfTopBranch = originalTree.branchLength(topology.getChildren().get(0).getContents(), topology.getChildren().get(1).getContents()) / 2.0;
      List<Taxon> leavesNames = new ArrayList<Taxon>();
      List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
      for (Arbre<CoalescentNode> current : model.roots)
      {
        if (!current.isLeaf()) throw new RuntimeException();
        Taxon t = current.getContents().nodeIdentifier;
        if (taxaToConsider.contains(t))
        {
          leavesNames.add(t);
          leaves.add(current.getContents().likelihoodModelCache);
        }
      }
      if (additionalLeaf != null)
      {
        leavesNames.add(additionalLeaf);
        leaves.add(additionalLLC);
      }
      current = PartialCoalescentState.initialState(
          leaves, 
          leavesNames,
          model.getObservations(),
          false);
      list.add(current);
    }
    
    private void iterate()
    {
      // pick a root
      Set<Taxon> fringe = set();
      for (Arbre<CoalescentNode> node : current.roots)
        fringe.add(node.getContents().nodeIdentifier);
      double [] prs = new double[current.nRoots()];
      for (int i = 0; i < prs.length; i++)
        prs[i] = fringe.contains(findParentAndSibling(i).getSecond()) ? 1 : 0;
      NumUtils.normalize(prs);
      int sampledIndex = SampleUtils.sampleMultinomial(rand, prs);
      final Taxon sampledTaxon = current.roots.get(sampledIndex).getContents().nodeIdentifier;
      Pair<Taxon,Taxon> parentAndSibling = findParentAndSibling(sampledIndex);
      final Taxon parent = parentAndSibling.getFirst(),
                  sibling= parentAndSibling.getSecond();
      int otherIndex = -1;
      aLoop:for (int i = 0; i < prs.length ; i++)
        if (current.roots.get(i).getContents().nodeIdentifier.equals(sibling))
        {
          otherIndex = i;
          break aLoop;
        }
      // find branch lengths
      double bl1 = parent == Taxon.dummy ? halfTopBranch : originalTree.branchLength(sampledTaxon, parent),
             bl2 = parent == Taxon.dummy ? halfTopBranch : originalTree.branchLength(sibling, parent);
      current = current.coalesce(sampledIndex, otherIndex, 0.0, bl1, bl2, parent);
      list.add(current);
    }
    
    private Pair<Taxon,Taxon> findParentAndSibling(int i)
    {
      Taxon pickedRoot = current.roots.get(i).getContents().nodeIdentifier;
      if (!taxaToConsider.contains(pickedRoot))
        throw new RuntimeException();
      // find its pointer in the guide
      Arbre<Taxon> pickedNode = Arbre.findFirstNodeWithContents(topology, pickedRoot);
      Arbre<Taxon> parent = pickedNode.getParent();
      if (!taxaToConsider.contains(parent.getContents()))
        throw new RuntimeException();
      for (Arbre<Taxon> child : parent.getChildren())
        if (taxaToConsider.contains(child.getContents()) && !child.getContents().equals(pickedNode.getContents()))
          return Pair.makePair(parent.getContents(), child.getContents());
      throw new RuntimeException();
    }
    
  }
  
  public static void main(String [] args)
  {
    NCPriorPriorKernel.deltaProposalRate = 10;
    //UnrootedTree ut = UnrootedTree.fromNewickRemovingBinaryRoot(new File("/Users/bouchard/w/legacy/state/remote/gp1187.seg1197.time1308586425516.exec/output/sim--1731337436.newick"));
    UnrootedTree ut = UnrootedTree.fromNewickRemovingBinaryRoot(new File("/Users/liangliang/experimental-results/bigtree/state/simpletree.newick"));
    System.out.println(ut);
    Taxon l1 = new Taxon("leaf_A"), l2 = new Taxon("leaf_E");
    CondSMCDissector cd = new CondSMCDissector(ut, l1, l2);
    Random rand = new Random(1);
    cd.sampleAnchor(rand);
    //
    PartialCoalescentState pcs = PartialCoalescentState.initState(ut.leaves(),false);
    List<PartialCoalescentState> list = cd.currentSampledPath(pcs, rand);
    
    StoreProcessor<PartialCoalescentState> pro = new StoreProcessor<PartialCoalescentState>();
    
    PartialCoalescentState init = list.get(0);
    
    List<PartialCoalescentState> path = list.subList(1, list.size());
    
    // set init state
    NCPriorPriorKernel kernel = new NCPriorPriorKernel(init);
    
    // compute weights
    double [] weights  = new double[path.size()];
//    for (int i = 0; i < weights.length; i++)
//      weights[i] =  (NCPriorPriorKernel.useOptimal ? path.get(i).nonClockLogWeight() : path.get(i).logLikelihoodRatio() - Math.log(path.get(i).nNonTrivialRoots()));
    
    ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
    
    // set the conditioning and its weights
    pf.setConditional(path, weights);
    
    // do the sampling
    pf.sample(kernel, pro);
    
    // sample from the samples of the pf
    PartialCoalescentState sampled = pro.sample(rand);
    
//    for (PartialCoalescentState cur : list)
//    {
//      System.out.println(cur);
//      System.out.println("====================");
//    }
//    System.out.println("Final );
//    System.out.println(RootedTree.Util.toString(list.get(list.size() -1).getFullCoalescentState()));
    
    System.out.println(cd.reconstitute(sampled));
//    System.out.println("anchor=" + cd.anchor);
//    System.out.println("san=" + cd.san);
  }
  
  
}
