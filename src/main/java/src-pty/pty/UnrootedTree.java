package pty;
import static nuts.util.CollUtils.map;
import fig.basic.UnorderedPair;
import goblin.DataPrepUtils;
import goblin.Taxon;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import ma.newick.NewickParser;
import ma.newick.ParseException;
import nuts.io.IO;
import nuts.math.Graph;
import nuts.math.Graphs;
import nuts.math.HashGraph;
import nuts.math.SemiGraph;
import nuts.util.Arbre;
import nuts.util.Arbre.ArbreMap;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.Indexer;
import nuts.util.MathUtils;
import nuts.util.Tree;
import pty.RootedTree.RootingInfo;
import pty.eval.SymmetricDiff;
import conifer.Phylogeny;
//import org.apache.lucene.analysis.ReusableAnalyzerBase;

public class UnrootedTree implements Phylogeny
{
  private static final long serialVersionUID = 1L;


  private final Graph<Taxon> 
    topo;
  

  private final Map<UnorderedPair<Taxon,Taxon>,Double> 
    branchLengths;
  
  /**
   * Note: use is not recommended (could behave badly with using Phylogeny objects)
   * @param edge
   * @param newValue
   */
  public void changeBranchLength(UnorderedPair<Taxon,Taxon> edge, double newValue)
  {
    if (!branchLengths.containsKey(edge))
      throw new RuntimeException();
    branchLengths.put(edge, newValue);
  }
  
//  private List<UnorderedPair<Taxon,Taxon>> 
//    edges = null,
//    nonTerminalEdge = null;
  
  public Graph<Taxon> getTopology() { return topo; }
  
  public List<Taxon> leaves()
  {
    List<Taxon> result = new ArrayList<Taxon>();
    for (Taxon l : topo.vertexSet())
      if (topo.nbrs(l).size() == 1)
        result.add(l);
    return result;
  }
  
  public Set<Taxon> leavesSet()
  {
    return new HashSet<Taxon>(leaves());
  }
  
  public UnrootedTree(UnrootedTree model) 
  {
    this.topo = model.topo;
    this.branchLengths = model.branchLengths;
  }
  
  public static UnrootedTree fromNewick(File path)
  {
    String treeDescr = IO.f2s(path);
    return fromNewick(treeDescr);
  }
  public static UnrootedTree fromNewickRemovingBinaryRoot(File f)
  {
    RootedTree t = RootedTree.Util.load(f);
    return UnrootedTree.fromRooted(t);
  }
  public static UnrootedTree fromNewickRemovingBinaryRoot(String stringContainingOneNewickFormattedTree)
  {
    RootedTree t = RootedTree.Util.fromNewickString(stringContainingOneNewickFormattedTree);
    return UnrootedTree.fromRooted(t);
  }
  public static UnrootedTree fromNewick(String stringContainingOneNewickFormattedTree)
  {
    try
    {
      final NewickParser np = new NewickParser(stringContainingOneNewickFormattedTree);
      final Arbre<String> t = Arbre.tree2Arbre(np.parse());
      final Map<Taxon,Double> bl = np.getBranchLengths();
      Map<UnorderedPair<Taxon,Taxon>,Double> convertedBL 
       = new HashMap<UnorderedPair<Taxon,Taxon>,Double>();
      final Set<Taxon> vertices = new HashSet<Taxon>();
      final Set<UnorderedPair<Taxon,Taxon>> edges = new HashSet<UnorderedPair<Taxon,Taxon>>();
      for (Arbre<String> node : t.nodes())
      {
        Taxon current = new Taxon(node.getContents());
        vertices.add(current);
        if (!node.isRoot())
        {
          Taxon parent = new Taxon(node.getParent().getContents());
          UnorderedPair<Taxon,Taxon> key = new UnorderedPair<Taxon,Taxon>(current,parent);
          edges.add(key);
          convertedBL.put(key, bl.get(current));
        }
      }
      return new UnrootedTree(new HashGraph<Taxon>(vertices,edges), convertedBL);
    }
    catch (Exception e) { throw new RuntimeException(e); }
  }
  
  /**
   * Also removes binary stems, if any
   * TODO: decouple this step
   */
  public static UnrootedTree fromRooted(RootedTree c)
  {
//    Graph<Language> topo = new Graph.HashGraph<Language>(Arbre.arbre2Tree(c.topology()));
    Map<UnorderedPair<Taxon,Taxon>,Double>  branchLengths
     = new HashMap<UnorderedPair<Taxon,Taxon>,Double>();
    Set<Taxon> languages = new HashSet<Taxon>(c.topology().nodeContents());
    
    Set<UnorderedPair<Taxon,Taxon>> edges = 
      new HashSet<UnorderedPair<Taxon,Taxon>>();
    for (Arbre<Taxon> node : c.topology().nodes())
      if (!node.isRoot())
      {
        UnorderedPair<Taxon,Taxon> edge = null;
        double bl = -1;
        if (node.getParent().isRoot() && node.getParent().getChildren().size() == 2)
        { // special case if we are at the root: connect that two children, marginalizing root node
          languages.remove(c.topology().getContents());
          List<Arbre<Taxon>> parentChildren = 
            node.getParent().getChildren();
          edge = 
            new UnorderedPair<Taxon,Taxon>(
                parentChildren.get(0).getContents(), 
                parentChildren.get(1).getContents());
          bl = 
            c.branchLengths().get(parentChildren.get(0).getContents()) +
            c.branchLengths().get(parentChildren.get(1).getContents());
        }
        else
        {
          edge = 
            new UnorderedPair<Taxon,Taxon>(
                node.getContents(), 
                node.getParent().getContents());
          bl = c.branchLengths().get(node.getContents());
        }
        edges.add(edge);
        branchLengths.put(edge, bl);
      }
    return new UnrootedTree(new HashGraph<Taxon>(languages, edges), branchLengths);
  }
  
  @Deprecated
  /**
   * Use unRootedClades instead!
   */
  public Set<Set<Taxon>> clades()
  {
    // find an internal node
    Taxon internal = null;
    loop:for (Taxon lang : topo.vertexSet())
      if (topo.nbrs(lang).size() > 1)
      {
        internal = lang;
        break loop;
      }
    if (internal == null) throw new RuntimeException();
    return SymmetricDiff.cladesFromUnrooted(Arbre.tree2Arbre(Graphs.toTree(topo, internal)));
  }
  
  public Set<Set<Taxon>> unRootedClades()
  {
    Set<Set<Taxon>> result = CollUtils.set();
    for (UnorderedPair<Set<Taxon>, Set<Taxon>> bipart : inducedBiPartitions2BranchMap().keySet())
    {
      result.add(bipart.getFirst());
      result.add(bipart.getSecond());
    }
    return result;
  }
  
  /** 
   * WARNING: THE COUNTS ARE BRANCH LENGTHS!!
   * @return
   */
  public Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> inducedBiPartitions2BranchMap()
  {
    Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> result = new Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>>();
    List<Taxon> allLeaves = leaves();
    Arbre<Taxon> topology = Arbre.tree2Arbre(Graphs.toTree(topo, allLeaves.get(0)));
    Map<Arbre<Taxon>, Set<Taxon>> leavesMap = Arbre.leavesMap(topology);
    for (Arbre<Taxon> key : leavesMap.keySet())
      if (!key.isRoot())
      {
        Set<Taxon> clade = leavesMap.get(key);
        double bl = branchLength(key.getContents(), key.getParent().getContents());
        Set<Taxon> complement = CollUtils.set(allLeaves);
        complement.removeAll(clade);
        result.setCount(new UnorderedPair<Set<Taxon>,Set<Taxon>>(complement, clade), bl);
      }
    return result;
  }
  public static double partitionMetric(UnrootedTree ut1, UnrootedTree ut2)
  {
    return partitionMetric(ut1, ut2, false, false);
  }  
  public static double normalizedPartitionMetric(UnrootedTree ut1, UnrootedTree ut2, boolean useTightNormalizer)
  {
    return partitionMetric(ut1, ut2, true, useTightNormalizer);
  }
  public static double normalizedPartitionMetric(UnrootedTree ut1, UnrootedTree ut2)
  {
    return partitionMetric(ut1, ut2, true, false);
  }
  public static double tightlyNormalizedPartitionMetric(UnrootedTree ut1, UnrootedTree ut2)
  {
    return partitionMetric(ut1, ut2, true, true);
  }

	public static double matchingMetric(UnrootedTree ref, UnrootedTree guess) {
		// TODO Auto-generated method stub
		return 0;
	}

  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  public static double partitionMetric(UnrootedTree ut1, UnrootedTree ut2, boolean normalize, boolean useTightNormalizer)
  {
    if (!ut1.leavesSet().equals(ut2.leavesSet()))
      throw new RuntimeException();
    Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>>
      biparts1 = ut1.inducedBiPartitions2BranchMap(),
      biparts2 = ut2.inducedBiPartitions2BranchMap();
    for (UnorderedPair<Set<Taxon>,Set<Taxon>> bipartition : CollUtils.union(biparts1.keySet(), biparts2.keySet()))
      if (normalize && useTightNormalizer &&
          (bipartition.getFirst().size() <= 1 || bipartition.getSecond().size() <= 1))
      {
        if (!biparts1.keySet().contains(bipartition) || !biparts2.keySet().contains(bipartition))
          throw new RuntimeException();
        // TIGHT NORMALIZATION:
        //ignore the trivial bipartitions: everybody have them!
        //this avoid having misleading normalized distances, e.g. if reporting 1-distance
        biparts1.setCount(bipartition, 0.0);
        biparts2.setCount(bipartition, 0.0);
      }
      else
      {
        if (biparts1.getCount(bipartition) != 0.0)
          biparts1.setCount(bipartition,1.0);
        if (biparts2.getCount(bipartition) != 0.0)
          biparts2.setCount(bipartition,1.0);
      }
    return bipartitionMetric(biparts1, biparts2, true, normalize);
  }
  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  public static double robinsonFouldsMetric(UnrootedTree ut1, UnrootedTree ut2)
  {
    return robinsonFouldsMetric(ut1, ut2, false);
  }
  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  public static double normalizedRobinsonFouldsMetric(UnrootedTree ut1, UnrootedTree ut2)
  {
    return robinsonFouldsMetric(ut1, ut2, true);
  }
  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  public static double robinsonFouldsMetric(UnrootedTree ut1, UnrootedTree ut2, boolean normalize)
  {
    if (!ut1.leavesSet().equals(ut2.leavesSet()))
      throw new RuntimeException();
    Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>>
      biparts1 = ut1.inducedBiPartitions2BranchMap(),
      biparts2 = ut2.inducedBiPartitions2BranchMap();
    return bipartitionMetric(biparts1, biparts2, true, normalize);
  }
  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  public static double kuhnerFelsenstein(UnrootedTree ut1, UnrootedTree ut2)
  {
    return kuhnerFelsenstein(ut1, ut2, false);
  }
  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  public static double normalizedKuhnerFelsenstein(UnrootedTree ut1, UnrootedTree ut2)
  {
    return kuhnerFelsenstein(ut1, ut2, true);
  }
  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  public static double kuhnerFelsenstein(UnrootedTree ut1, UnrootedTree ut2, boolean normalize)
  {
    if (!ut1.leavesSet().equals(ut2.leavesSet()))
      throw new RuntimeException();
    Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>>
      biparts1 = ut1.inducedBiPartitions2BranchMap(),
      biparts2 = ut2.inducedBiPartitions2BranchMap();
    return bipartitionMetric(biparts1, biparts2, false, normalize);
  }
  /**
   * Warning: there might be slight numerical difference across calls depending on the order the bipartitions set are enumerated
   */
  private static double bipartitionMetric(
      Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> biparts1, 
      Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> biparts2, boolean RF, boolean normalize)
  {
    double n1 = 1.0, n2 = 1.0;
    if (normalize)
    {
      n1 = biparts1.totalCount() * 2;
      n2 = biparts2.totalCount() * 2;
    }
    double sum = 0.0;
    for (UnorderedPair<Set<Taxon>,Set<Taxon>> bipartition : CollUtils.union(biparts1.keySet(), biparts2.keySet()))
    {
      final double b1 = biparts1.getCount(bipartition) / n1,
                   b2 = biparts2.getCount(bipartition) / n2;
      final double diff = b1 - b2;
      sum += Math.abs( diff * (RF ? 1.0 : diff));
    }
    return sum;
  }
  
  public double diameter()
  {
	  Counter<UnorderedPair<Taxon,Taxon>> pairwiseEdgeNumbers=pairwiseEdgeNumberDistances();
//	  System.out.println(pairwiseEdgeNumbers);
	  return pairwiseEdgeNumbers.max();
  }
  
  public Counter<UnorderedPair<Taxon,Taxon>> pairwiseEdgeNumberDistances()
  {
    Counter<UnorderedPair<Taxon,Taxon>> result = new Counter<UnorderedPair<Taxon,Taxon>>();
    List<Taxon> leaves = leaves();    
    for (int i = 0; i < leaves.size(); i++)
      for (int j = i+1; j < leaves.size(); j++)
      {
        Taxon l1 = leaves.get(i), l2 = leaves.get(j);
        result.setCount(new UnorderedPair<Taxon,Taxon>(l1,l2), pairwiseEdgeNumberDistance(l1, null, l2));
      }
    return result;
  }
  
  private double pairwiseEdgeNumberDistance(Taxon current, Taxon parent, Taxon target)
  {
    if (current.equals(target)) return 0;
    double min = Double.POSITIVE_INFINITY;
    for (Taxon nb : topo.nbrs(current))
      if (parent == null || !nb.equals(parent))
      {
        final double recur = pairwiseEdgeNumberDistance(nb, current, target);
        if (!Double.isInfinite(recur)) 
          return 1 + recur;
      }
    return min;
  }
  
   
  public Counter<UnorderedPair<Taxon,Taxon>> pairwiseDistances()
  {
    Counter<UnorderedPair<Taxon,Taxon>> result = new Counter<UnorderedPair<Taxon,Taxon>>();
    List<Taxon> leaves = leaves();    
    for (int i = 0; i < leaves.size(); i++)
      for (int j = i+1; j < leaves.size(); j++)
      {
        Taxon l1 = leaves.get(i), l2 = leaves.get(j);
        result.setCount(new UnorderedPair<Taxon,Taxon>(l1,l2), pairwiseDistance(l1, null, l2));
      }
    return result;
  }
  
  private double pairwiseDistance(Taxon current, Taxon parent, Taxon target)
  {
    if (current.equals(target)) return 0.0;
    double min = Double.POSITIVE_INFINITY;
    for (Taxon nb : topo.nbrs(current))
      if (parent == null || !nb.equals(parent))
      {
        final double recur = pairwiseDistance(nb, current, target);
        if (!Double.isInfinite(recur)) 
          return branchLength(nb, current) + recur;
      }
    return min;
  }
  public Tree<Taxon> toTree(Taxon aNode)
  {
    return Graphs.toTree(topo, aNode);
  }
  public Tree<Taxon> toTree()
  {
    return Graphs.toTree(topo, leaves().get(0));
  }
  /**
   * return null if the rerooting is not valid (ie. l1-l2 is no longer an edge)
   * @param rf
   * @return
   */
  public RootedTree reRoot(RootingInfo rf)
  {
    return reRoot(rf.root, rf.l1, rf.l2, rf.ratioToL1);
  }
  private RootedTree reRoot(Taxon rootName, Taxon l1, Taxon l2, double ratioToL1)
  {
    if (!MathUtils.isCloseToProb(ratioToL1))
      throw new RuntimeException();
    Tree<Taxon> tree = toTree(l1);
    Tree<Taxon> subt1 = new Tree<Taxon>(l1),
                   subt2 = null;
    for (Tree<Taxon> st : tree.getChildren())
      if (st.getLabel().equals(l2))
        subt2 = st;
      else
        subt1.getChildren().add(st);
    
    Tree<Taxon> fullTree = new Tree<Taxon>(rootName);
    fullTree.getChildren().add(subt1);
    if (subt2 == null) 
      return null;
    fullTree.getChildren().add(subt2);
    Arbre<Taxon> a = Arbre.tree2Arbre(fullTree);
    // now, branch lengths
    final Map<Taxon,Double> bls = CollUtils.map();
    for (Arbre<Taxon> sa : a.nodes())
      if (sa.isRoot()) 
        ;
      else 
      {
        final Taxon cur = sa.getContents();
        double val = Double.NaN;
        if (cur.equals(l1))
          val = ratioToL1 * branchLength(l1,l2);
        else if (cur.equals(l2))
          val = (1.0 - ratioToL1) * branchLength(l1,l2);
        else
          val = branchLength(cur, sa.getParent().getContents());
        bls.put(cur, val);
      }
    return RootedTree.Util.create(a, bls);
  }
  
//  public void restrict(Set<Language> toRetain)
//  {
//    Language aLang = CollUtils.pick(topo.vertexSet());
//    Language aNei  = CollUtils.pick(topo.nbrs(aLang));
//    RootedTree rooted = reRoot(new RootingInfo(aLang, aNei, new Language("NewRoot"), 0.5));
//    Tree<Language> t = Arbre.arbre2Tree(rooted.topology());
//    t = restrict(t, toRetain);
//    t = removeUselessInternalNodes(t);
//    System.out.println(t);
//  }
  

  
  public String toNewick()
  {
    Tree<Taxon> t = Graphs.toTree(topo);
    return DataPrepUtils.newick2(t, branchLengths, false);
  }
  
  public String toNewick(Taxon root)
  {
    Tree<Taxon> t = Graphs.toTree(topo, root);
    return DataPrepUtils.newick2(t, branchLengths, false);
  }
  
  public List<UnorderedPair<Taxon,Taxon>> nbrEdges(UnorderedPair<Taxon,Taxon> edge)
  {
    List<UnorderedPair<Taxon,Taxon>> result = new ArrayList<UnorderedPair<Taxon,Taxon>>();
    for (Taxon lang : topo.nbrs(edge.getFirst ())) result.add(new UnorderedPair<Taxon,Taxon>(edge.getFirst() , lang));
    for (Taxon lang : topo.nbrs(edge.getSecond())) result.add(new UnorderedPair<Taxon,Taxon>(edge.getSecond(), lang));
    return result;
  }
  
  @Override
  public String toString()
  {
    Arbre<Taxon> a = Arbre.tree2Arbre(Graphs.toTree(topo));
    Arbre<String> annotated = a.preOrderMap(new ArbreMap<Taxon,String>() {
      @Override
      public String map(Arbre<Taxon> currentDomainNode)
      {
        String result = currentDomainNode.getContents().toString();
        if (!currentDomainNode.isRoot())
          result += ":" + branchLength(
              currentDomainNode.getContents(), 
              currentDomainNode.getParent().getContents());
        return result;
      }
    });
    return annotated.deepToString();
  }  
  
  public static UnrootedTree normalize(UnrootedTree ut)
  {
    RootedTree rt = RootedTree.Util.centroidRooting(ut);
    rt = RootedTree.Util.normalizeBranches(rt);
    return fromRooted(rt);
  }
  
  public UnrootedTree(Graph<Taxon> topo, Map<UnorderedPair<Taxon,Taxon>,Double> branchLengths)
  {
    this.topo = topo;
    this.branchLengths = branchLengths;
  }
  
  public double branchLength(Taxon l1, Taxon l2) 
  {
    return branchLength(new UnorderedPair<Taxon,Taxon>(l1,l2));
  }
  public double branchLength(UnorderedPair<Taxon,Taxon> edge) 
  {
    return branchLengths.get(edge);
  }
  
  public static class EfficientUnrootedTree
  {
    private final int [][] nbhrs;
    private final double [][] bls;
    private final int nLeaves;
    private final Indexer<Taxon> indexer = new Indexer<Taxon>();
    
    public EfficientUnrootedTree(UnrootedTree ut)
    {
      final List<Taxon> leaves = ut.leaves();
      this.nLeaves = leaves.size();
      for (final Taxon leaf : leaves)
        indexer.addToIndex(leaf);
      for (final Taxon node : ut.getTopology().vertexSet())
        if (!indexer.containsObject(node))
          indexer.addToIndex(node);
      final int size = indexer.size();
      this.nbhrs = new int[size][];
      this.bls = new double[size][];
      for (final Taxon t1 : ut.getTopology().vertexSet())
      {
        final int i1 = indexer.o2i(t1);
        Set<Taxon> curNbhrs = ut.topo.nbrs(t1);
        nbhrs[i1] = new int[curNbhrs.size()];
        bls[i1] = new double[curNbhrs.size()];
        int idx = 0;
        for (Taxon t2 : ut.getTopology().nbrs(t1))
        {
          final int i2 = indexer.o2i(t2);
          nbhrs[i1][idx] = i2;
          bls[i1][idx] = ut.branchLength(t1, t2);
          idx++;
        }
      }
    }
    
    public Counter<UnorderedPair<Taxon,Taxon>> allTotalBranchLengthDistances()
    {
      double [][] result = new double[nLeaves][nLeaves];
      
      for (int start = 0; start < nLeaves - 1; start++)
        _dfsTotalBL(0.0, start, start, -1, result);
      
      // conversion
      Counter<UnorderedPair<Taxon,Taxon>> convertedResult = new Counter<UnorderedPair<Taxon,Taxon>>();
      for (int l1 = 0; l1 < nLeaves; l1++)
      {
        final Taxon t1 = indexer.i2o(l1);
        for (int l2 = l1 + 1; l2 < nLeaves; l2++)
          convertedResult.setCount(new UnorderedPair<Taxon, Taxon>(t1, indexer.i2o(l2)), result[l1][l2]);
      }
      return convertedResult;
    }

    private void _dfsTotalBL(double parentLen, int start, int current, int parent, double[][] result)
    {
      final int [] thisNbhrs = nbhrs[current];
      final double [] thisBLs = bls[current];
      if (thisNbhrs.length != thisBLs.length)
        throw new RuntimeException();
      for (int nIndex = 0; nIndex < thisNbhrs.length; nIndex++)
      {
        int nbhr = thisNbhrs[nIndex];
        if (nbhr != parent)
        {
          final double newLen = parentLen + thisBLs[nIndex];
          if (nbhr < nLeaves)
            // a leaf!
            result[start][nbhr] = newLen;
          else
            // recurse!
            _dfsTotalBL(newLen, start, nbhr, current, result);
        }
      }
    }
  }
  
  public Counter<UnorderedPair<Taxon,Taxon>> allTotalBranchLengthDistances()
  {
    return new EfficientUnrootedTree(this).allTotalBranchLengthDistances();
    
//    Counter<UnorderedPair<Taxon,Taxon>> result = new Counter<UnorderedPair<Taxon,Taxon>>();
//    Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> bm = inducedBiPartitions2BranchMap();
//    LogInfo.track("Method 1");
//    for (Entry<UnorderedPair<Set<Taxon>, Set<Taxon>>, Double> bipartPair : bm.entries.entrySet())
////    for (UnorderedPair<Set<Taxon>, Set<Taxon>> bipart : bm.keySet())
//    {
//      UnorderedPair<Set<Taxon>, Set<Taxon>> bipart = bipartPair.getKey();
//      double count = bipartPair.getValue();
//      for (Taxon t1 : bipart.getFirst())
//        for (Taxon t2 : bipart.getSecond())
//          result.incrementCount(new UnorderedPair<Taxon, Taxon>(t1, t2), count);
//    }
//    LogInfo.end_track();
//    
//    LogInfo.track("Method 2");
//    
//    Counter<UnorderedPair<Taxon,Taxon>> result2 = new EfficientUnrootedTree(this).allTotalBranchLengthDistances();
//    
//    LogInfo.end_track();
//    
//    for (UnorderedPair<Taxon,Taxon> key : CollUtils.union(result.keySet(), result2.keySet()))
//      MathUtils.checkClose(result.getCount(key), result2.getCount(key));
//    
//    
////    System.out.println("CHecking..");
////    Counter test = new Counter();
////    List<Taxon> taxa = leaves();
////  for (int i = 0 ; i < taxa.size(); i++)
////  for (int j = i+1; j< taxa.size() ;j ++)
////  {
////    Taxon t1 = taxa.get(i), t2 = taxa.get(j);
////    test.incrementCount(new UnorderedPair<Taxon,Taxon>(t1,t2), this.totalBranchLengthDistance(t1,t2) * w);
////  }
////  for (Object key : CollUtils.union(test.keySet(), result.keySet()))
////  {
////    double v1 = test.getCount(key), v2 = result.getCount((UnorderedPair) key);
////    if (!MathUtils.close(v1, v2))
////      System.out.println("" + v1 + " vs " + v2);
////  }
//  
//  
//    return result;
  }
  
  public double totalBranchLength()
  {
    double sum = 0.0;
    for (UnorderedPair<Taxon, Taxon> edge : edges())
      sum += branchLength(edge);
    return sum;
  }
  
  public double totalBranchLengthDistance(Taxon l1, Taxon l2)
  {
    Tree<Taxon> t = toTree(l1);
    return _totalBranchLengthDist(t, l2);
  }
  
  private double _totalBranchLengthDist(Tree<Taxon> t, Taxon target)
  {
    if (t.getLabel().equals(target)) return 0.0;
    double result = Double.POSITIVE_INFINITY;
    for (Tree<Taxon> child : t.getChildren())
    {
      double cur = _totalBranchLengthDist(child, target);
      if (!Double.isInfinite(cur))
      {
        if (!Double.isInfinite(result)) throw new RuntimeException();
        result = cur + branchLength(t.getLabel(), child.getLabel());
      }
    }
    return result;
  }

  public UnorderedPair<Taxon,Taxon> randomEdge(Random rand)
  {
    return edges().get(rand.nextInt(edges().size()));
  }
  
  public UnorderedPair<Taxon,Taxon> randomNonTerminalEdge(Random rand)
  {
//    if (nonTerminalEdges().size() <= 0)
//    {
//      System.out.println("Edges: " + edges());
//      for (UnorderedPair<Language,Language> edge : edges())
//        System.out.println ("\tEdge:" + edge + " nbrs: " + topo.nbrs(edge.getFirst()) + " and: " + topo.nbrs(edge.getSecond()));
//    }
    if (nonTerminalEdges().size() == 0)
      return null;
    return nonTerminalEdges().get(rand.nextInt(nonTerminalEdges().size()));
  }
  
  public List<UnorderedPair<Taxon, Taxon>> edges()
  {
    return _deterministic_edges();
//    if (edges == null)
//      edges = _deterministic_edges(); //new ArrayList<UnorderedPair<Taxon,Taxon>>(Graphs.edgeSet(topo));
//    return edges;
  }
  
  private List<UnorderedPair<Taxon, Taxon>> _deterministic_edges()
  {
    
    List<UnorderedPair<Taxon,Taxon>> result = CollUtils.list();
    List<Taxon> taxa = CollUtils.list(topo.vertexSet());
    Collections.sort(taxa);
    for (Taxon t : taxa)
    {
      List<Taxon> nbs = CollUtils.list(topo.nbrs(t));
      Collections.sort(nbs);
      for (Taxon other : nbs)
        if (other.compareTo(t) > 0)
          result.add(new UnorderedPair<Taxon,Taxon>(t, other));
    }
    
//    {
//      IO.warnOnce("EXPENSIVE CHECK BEING DONE");
//      
//      if (!CollUtils.set(result).equals(CollUtils.set(_non_deterministic_edges())))
//        throw new RuntimeException();
//    }
    
    return result;
  }
  
  private List<UnorderedPair<Taxon, Taxon>> _non_deterministic_edges()
  {
    return new ArrayList<UnorderedPair<Taxon,Taxon>>(Graphs.edgeSet(topo));
  }
  
  public List<UnorderedPair<Taxon, Taxon>> nonTerminalEdges()
  {
    List<UnorderedPair<Taxon, Taxon>> nonTerminalEdge = null;
    if (nonTerminalEdge == null)
    {
      nonTerminalEdge = new ArrayList<UnorderedPair<Taxon,Taxon>>();
      for (UnorderedPair<Taxon,Taxon> edge : edges())
        if (topo.nbrs(edge.getFirst()).size() > 1 &&
            topo.nbrs(edge.getSecond()).size() > 1)
          nonTerminalEdge.add(edge);
    }
    return nonTerminalEdge;
  }
  
  public static UnrootedTree removeZeroes(UnrootedTree ut)
  {
    Map<UnorderedPair<Taxon, Taxon>,Double> bls = map();
    double minNotZero = Double.POSITIVE_INFINITY;
    for (UnorderedPair<Taxon, Taxon> bl : ut.branchLengths.keySet())
    {
      final double curLen = ut.branchLengths.get(bl);
      if (curLen> 0 && curLen < minNotZero)
        minNotZero = curLen;
    }
    for (UnorderedPair<Taxon, Taxon> bl : ut.branchLengths.keySet())
    {
      final double curLen = ut.branchLengths.get(bl);
      if (curLen == 0.0)
        bls.put(bl, minNotZero/2.0);
      else
        bls.put(bl, curLen);
    }
    return new UnrootedTree(ut.topo, bls);
  }
  
  
  // copy and modify for mcmc moves
  
  public UnrootedTree branchLengthNeighbor(UnorderedPair<Taxon, Taxon> edge, double newLength)
  {
    Map<UnorderedPair<Taxon,Taxon>,Double> branchLengthsCopy =
      new HashMap<UnorderedPair<Taxon,Taxon>,Double>(branchLengths);
    branchLengthsCopy.put(edge, newLength);
    return new UnrootedTree(this.topo, branchLengthsCopy);
  }
  
  public List<UnrootedTree> topologicalNeighbors(UnorderedPair<Taxon, Taxon> edge)
  {
    List<UnrootedTree> result = new ArrayList<UnrootedTree>();
    // name the endpoints arbitrarily
    final Taxon 
      x = edge.getFirst(),
      v = edge.getSecond();
    // get the nbrs excluding x,v:
    final List<Taxon>
      xN = new ArrayList<Taxon>(),
      vN = new ArrayList<Taxon>();
    for (Taxon l : topo.nbrs(x)) if (l != v) xN.add(l);
    for (Taxon l : topo.nbrs(v)) if (l != x) vN.add(l);
    if (xN.size() != 2 || vN.size() != 2) throw new RuntimeException();
    // name one of the other neighbors of v arbitrarily (that is, other than x)
    final Taxon y = vN.get(0);
    Collections.sort(xN);
    for (Taxon z : xN)
      result.add(neighborInterchange(v,x,y,z));
    return result;
  }
  /**
   * In Newick notation:
   * take the tree (((...)z,(...)a)x,((...)y,(...)b)v); 
   * and make it   (((...)y,(...)a)x,((...)z,(...)b)v);
   * @param v
   * @param x
   * @param y
   * @param z
   * @return
   */
  private UnrootedTree neighborInterchange(
      final Taxon v, 
      final Taxon x, 
      final Taxon y,
      final Taxon z)
  {
    final UnorderedPair<Taxon,Taxon>
      zx = new UnorderedPair<Taxon,Taxon>(z, x),
      yv = new UnorderedPair<Taxon,Taxon>(y, v),
      zv = new UnorderedPair<Taxon,Taxon>(z, v),
      yx = new UnorderedPair<Taxon,Taxon>(y, x);
    // create the graph
    Graph<Taxon> newTopo = new HashGraph<Taxon>(new SemiGraph<Taxon>() {
      public boolean hasSemiEdge(Taxon one, Taxon two)
      {
        final UnorderedPair<Taxon,Taxon> cur = new UnorderedPair<Taxon,Taxon>(one, two);
             if (cur.equals(zx) || cur.equals(yv)) return false;
        else if (cur.equals(zv) || cur.equals(yx)) return true;
        else return topo.hasEdge(one, two);
      }
      public Set<Taxon> vertexSet() { return topo.vertexSet(); }
    });
    // create the new branch lengths
    Map<UnorderedPair<Taxon,Taxon>,Double> newBranches = 
      new HashMap<UnorderedPair<Taxon,Taxon>,Double>();
    newBranches.put(yx, branchLengths.get(yv));
    newBranches.put(zv, branchLengths.get(zx));
    for (UnorderedPair<Taxon,Taxon> branch : branchLengths.keySet())
      if (!branch.equals(zx) && !branch.equals(yv))
        newBranches.put(branch, branchLengths.get(branch));
    return new UnrootedTree(newTopo, newBranches);
  }
  
  public static void main(String [] args) throws ParseException, IOException
  {
    {
      UnrootedTree temp = UnrootedTree.fromNewick(new File("e/575.exec/WALS-chain-0/mle.newick"));
      System.out.println("Original unrooted:\n" + temp);
      RootingInfo rt = new RootingInfo(new Taxon("internal_4"), new Taxon("internal_5"), new Taxon("root"), 0.5);
      RootedTree t = temp.reRoot(rt);
      System.out.println("Rerooted:\n" + t);
      RootingInfo t2 = RootedTree.Util.getRootingInfo(t);
      UnrootedTree temp2 = UnrootedTree.fromRooted(t);
      System.out.println("Rerooted again:\n" + temp2.reRoot(t2) );
    
//      RootedTree rt = RootedTree.Util.load(new File("e/575.exec/WALS-chain-0/mle.newick"));
//      System.out.println("Original rooted:\n" + rt);
//      RootingInfo ri = RootedTree.Util.getRootingInfo(rt);
//      NonClockTree nct = NonClockTree.fromCoalescent(rt);
//      System.out.println("Transformed into unrooted:\n" + nct);
//      System.out.println("Rerooted:\n" + nct.reRoot(ri));
    
    
    }
    
    
    // test reader
    UnrootedTree temp = UnrootedTree.fromNewick(new File("e/575.exec/WALS-chain-0/mle.newick"));
    System.out.println(temp);
    
    
    
    
    System.out.println(temp.leaves().size());
    System.out.println("Size:" +temp.clades().size());
    for (Set<Taxon> clade : temp.clades())
      System.out.println(clade);
    
    // get a coalescent
    RootedTree c = RootedTree.Util.load(new File("data/toy-tree.newick"));
    // convert into a nonclock representation
    UnrootedTree nct = UnrootedTree.fromRooted(c);
    System.out.println(nct);
    // do a few branch length resizing
    Random rand = new Random(1);
    for (int i = 0; i < 10000; i++)
    {
      UnorderedPair<Taxon,Taxon> edge = nct.randomEdge(rand);
      int newL = rand.nextInt(10);
//      System.out.println("Edge resampled: " + edge + ", new length: " + newL);
      nct = nct.branchLengthNeighbor(edge, newL);
//      System.out.println(nct);
      // do a few neighbor interchanges
      edge = nct.randomNonTerminalEdge(rand);
//      System.out.println("Edge to do stNNI: " + edge);
      for (UnrootedTree t : nct.topologicalNeighbors(edge))
      {
        nct = t;
//        System.out.println(t);
      }
    }
    System.out.println("done");
  }
  
  public static interface UnrootedTreeProcessor
  {
    public void process(UnrootedTree ut);
  }

  @Override
  public int nTaxa()
  {
    return leaves().size();
  }

  @Override
  public RootedTree getRooted()
  {
    return null;
  }

  @Override
  public UnrootedTree getUnrooted()
  {
    return this;
  }

}
