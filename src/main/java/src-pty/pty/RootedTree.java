package pty;
import static nuts.util.CollUtils.map;
import fig.basic.IOUtils;
import fig.basic.UnorderedPair;
import goblin.CognateId;
import goblin.DataPrepUtils;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import conifer.Phylogeny;

import ma.BalibaseCorpus;
import ma.BioCorpus;
import ma.newick.NewickParser;
import ma.newick.ParseException;
import ma.newick.TreeNode;
import nuts.math.Sampling;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Tree;
import nuts.util.Arbre.ArbreMap;
import pty.smc.models.CTMC;

public interface RootedTree extends Phylogeny
{
  public Arbre<Taxon> topology();
  public Map<Taxon,Double> branchLengths(); 
  
  public static class RootingInfo
  {
    public final Taxon l1, l2, root;
    public final double ratioToL1;
    public RootingInfo(Taxon l1, Taxon l2, Taxon root, double ratioToL1)
    {
      if (ratioToL1 < 0 || ratioToL1 > 1)
        throw new RuntimeException();
      this.l1 = l1;
      this.l2 = l2;
      this.root = root;
      this.ratioToL1 = ratioToL1;
    }
    @Override
    public String toString()
    {
      return "(" + root + " (" + l1 + ":" + ratioToL1 + "x ...) (" + l2 + ":" + (1.0 - ratioToL1) + "x ...))";
    }
  }

  public static interface RootedTreeProcessor
  {
    public void process(RootedTree rt);
  }
  
  public static class Util
  {
    
    public static class RootedTreeImpl implements RootedTree
    {
      private static final long serialVersionUID = 1L;
      private final Arbre<Taxon> topo;
      private final Map<Taxon,Double> bl;
      public RootedTreeImpl(Arbre<Taxon> topo, Map<Taxon, Double> bl)
      {
        this.topo = topo;
        this.bl = bl;
        if (bl.keySet().contains(topo.root().getContents()))
          throw new RuntimeException();
        Set<Taxon> nodes = CollUtils.set( bl.keySet() );
        nodes.add(topo.getContents());
        if (!nodes.equals(CollUtils.set(topo.nodeContents())))
          throw new RuntimeException();
      }
      @Override
      public Map<Taxon, Double> branchLengths() { return Collections.unmodifiableMap(bl); }
      @Override
      public Arbre<Taxon> topology() { return topo; } 
      @Override
      public String toString()
      {
        return Util.toString(this);
      }
      @Override
      public int nTaxa()
      {
        return topo.nLeaves();
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
    }
    
    public static RootedTree translate(RootedTree rt, Map<Taxon,Taxon> translation)
    {
      Arbre<Taxon> translatedA = Arbre.map(rt.topology(), translation);
      Map<Taxon,Double> translatedBL = map();
      for (Taxon t : rt.branchLengths().keySet())
        translatedBL.put(translation.get(t), rt.branchLengths().get(t));
      return new RootedTreeImpl(translatedA, translatedBL);
    }
    
    public static RootedTree normalizeBranches(RootedTree rt)
    {
      double sum = 0.0;
      for (Taxon t : rt.branchLengths().keySet())
        sum += rt.branchLengths().get(t);
      Map<Taxon,Double> newBLs = CollUtils.map();
      for (Taxon t : rt.branchLengths().keySet())
        newBLs. put(t, rt.branchLengths().get(t) / sum);
      return new RootedTreeImpl(rt.topology(), newBLs);
    }
    
    public static RootedTree coalesce(Taxon newRoot, RootedTree rt1, RootedTree rt2, double bl1, double bl2)
    {
      Arbre<Taxon> 
        left = rt1.topology().copy(),
        right= rt2.topology().copy();
      Arbre<Taxon> newTree = Arbre.arbreWithChildren(newRoot, left, right);
      Map<Taxon,Double> bls = CollUtils.map();
      bls.putAll(rt1.branchLengths());
      bls.putAll(rt2.branchLengths());
      bls.put(left.getContents(), bl1);
      bls.put(right.getContents(), bl2);
      return new RootedTreeImpl(newTree, bls);
    }
    
    public static RootedTreeImpl incrementSmallBranches(RootedTree rt, double min)
    {
      Map<Taxon,Double> incrementSmallBranches =  incrementSmallBranches(rt.branchLengths(), min);
      return new RootedTreeImpl(rt.topology(), incrementSmallBranches);
    }
    
    public static Map<Taxon,Double> incrementSmallBranches(Map<Taxon,Double> bls, double min)
    {
      Map<Taxon,Double> result = CollUtils.map();
      for (Taxon l : bls.keySet())
        if (bls.get(l) < min) result.put(l,min);
        else                  result.put(l,bls.get(l));
      return result;
    }
    
    public static RootedTree restrict(RootedTree rt, Set<Taxon> toRetain)
    {
      Tree<Taxon> t = Arbre.arbre2Tree(rt.topology());
      t = restrict(t, toRetain);
      t = removeUselessInternalNodes(t);
      // compute new branch lengths
      Map<Taxon,Double> bls = CollUtils.map();
      UnrootedTree originalUT = UnrootedTree.fromRooted(rt);
      Arbre<Taxon> a = Arbre.tree2Arbre(t);
      for (Arbre<Taxon> subt : a.nodes())
        if (!subt.isRoot())
          bls.put(subt.getContents(), 
              originalUT.totalBranchLengthDistance(subt.getContents(), subt.getParent().getContents()));
      return new RootedTreeImpl(a, bls);
    }
    
    public static Tree<Taxon> restrict(Tree<Taxon> t, Set<Taxon> toRetain)
    {
      if (t.isLeaf()) return (toRetain.contains(t.getLabel()) ? t : null);
      List<Tree<Taxon>> newChildren = new ArrayList<Tree<Taxon>>();
      Tree<Taxon> cur = null;
      for (Tree<Taxon> child : t.getChildren())
        if ((cur = restrict(child, toRetain)) != null)
            newChildren.add(cur);
      return newChildren.size() == 0 ? null : new Tree<Taxon>(t.getLabel(), newChildren);
    }
    public static Tree<Taxon> removeUselessInternalNodes(Tree<Taxon> t)
    {
      List<Tree<Taxon>> newChildren = new ArrayList<Tree<Taxon>>();
      if (t.getChildren().size() == 1) return removeUselessInternalNodes(t.getChildren().get(0));
      else
        for (Tree<Taxon> child : t.getChildren())
          newChildren.add(removeUselessInternalNodes(child));
      return new Tree<Taxon>(t.getLabel(), newChildren);
    }
    
    public static RootingInfo getRootingInfo(RootedTree rt)
    {
      if (rt.topology().getChildren().size() != 2)
        throw new RuntimeException();
      final Taxon 
        l1 = rt.topology().getChildren().get(0).getContents(),
        l2 = rt.topology().getChildren().get(1).getContents();
      final double b1 = rt.branchLengths().get(l1),
                   b2 = rt.branchLengths().get(l2);
      return new RootingInfo(
          l1,
          l2,
          rt.topology().getContents(),
          b1 / (b1 + b2));
    }
    
    public static RootedTree create(Arbre<Taxon> topo, Map<Taxon, Double> bl)
    {
      return new RootedTreeImpl(topo,bl);
    }
    
    public static String toString(final RootedTree nct)
    {
      return nct.topology().preOrderMap(new ArbreMap<Taxon,String>() {

        @Override
        public String map(Arbre<Taxon> currentDomainNode)
        {
          return currentDomainNode.getContents().toString() + 
            (currentDomainNode.isRoot() ? "" : ":" +  nct.branchLengths().get(currentDomainNode.getContents()));
        }
      }).deepToString();
    }
    
    public static String toNewick(RootedTree rt)
    {
      return DataPrepUtils.newick(Arbre.arbre2Tree(rt.topology()), rt.branchLengths(), false);
    }
    
    public static RootedTree fromBalibase(BioCorpus bc, CognateId id)
    {
      Tree<String> rawTopo = bc.getTopology(id);
      Arbre<String> converted = Arbre.tree2Arbre(rawTopo);
      Arbre<Taxon> converted2 = Taxon.LanguageUtils.convert(converted);
      return new RootedTreeImpl(converted2, bc.getBranchLengths(id));
    }
    
    public static RootedTree load(File pathToNewick)
    {
      BufferedReader reader = IOUtils.openInHard(pathToNewick);
      NewickParser np = new NewickParser(reader);
      try { 
        RootedTree result =  new RootedTreeImpl(Taxon.LanguageUtils.convert(Arbre.tree2Arbre(np.parse())), np.getBranchLengths());
        return result;
      } catch (ParseException e) { throw new RuntimeException(e); }
      finally { 
        try { reader.close(); } 
        catch (IOException e) { throw new RuntimeException(); } 
      }
    }
    public static RootedTree fromNewickString(String newickStr)
    {
      NewickParser np = new NewickParser(newickStr);
      
      try {
    	  return new RootedTreeImpl(Taxon.LanguageUtils.convert(Arbre.tree2Arbre(np.parse())), np.getBranchLengths());
      } catch (ParseException e) { throw new RuntimeException(e); }
    }
    /**
     * assumes an ultrametric tree
     */
    public static double height(RootedTree rt)
    {
      double sum = 0;
      Arbre<Taxon> cur = rt.topology().leaves().iterator().next();
      while (!cur.isRoot())
      {
        sum += rt.branchLengths().get(cur.getContents());
        cur = cur.getParent();
      }
      return sum;
    }

    /**
     * Weird random model.
     * Creates rooted tree with each branch exp(1)
     * 
     * Note: not the same as unrooted with all branch exp(1) b/c of rooting has 2 branches of exp(1)
     */
    public static RootedTree random(Random rand, Collection<Taxon> leaves)
    {
      Set<RootedTree> toMerge = CollUtils.set();
      for (Taxon leaf : leaves)
        toMerge.add(singleton(leaf));
      int i = 0;
      while (toMerge.size() > 1)
      {
        List<RootedTree> pair = CollUtils.list(Sampling.sampleSubset(rand, toMerge, 2));
        toMerge.removeAll(pair);
        final double length0 = Sampling.sampleExponential(rand, 1.0),
        length1 = Sampling.sampleExponential(rand, 1.0);
        final Taxon internal = new Taxon("internal_" + (i++));
        toMerge.add(coalesce(internal, pair.get(0), pair.get(1), length0, length1));
      }
      return CollUtils.pick(toMerge);
    }

    public static RootedTree singleton(Taxon t)
    {
      Arbre<Taxon> topo = Arbre.arbre(t);
      Map<Taxon,Double> bls = map();
      return RootedTree.Util.create(topo, bls);
    }

    public static RootedTree centroidRooting(UnrootedTree t)
    {
      // pick the edge
      UnorderedPair<Taxon,Taxon> edge = null;
      double min = Double.POSITIVE_INFINITY;
      for (UnorderedPair<Taxon,Taxon> e : t.edges())
      {
        double cur = sumHeights(t, e.getFirst()) + sumHeights(t, e.getSecond());
        if (cur < min)
        {
          edge = e;
          min = cur;
        }
      }
      RootingInfo rooting = new RootingInfo(edge.getFirst(), edge.getSecond(), new Taxon("root"),0.5);
      return t.reRoot(rooting);
    }

    private static double sumHeights(UnrootedTree t, Taxon first)
    {
      double sum = 0;
      for (Taxon leaf : t.leaves())
        sum += t.totalBranchLengthDistance(leaf, first);
      return sum;
    }
  }
}
