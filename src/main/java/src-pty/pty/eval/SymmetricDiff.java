package pty.eval;
import fig.basic.IOUtils;
import goblin.Taxon;
import goblin.BayesRiskMinimizer.LossFct;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;

import pty.smc.PartialCoalescentState;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.test.SymmetryTest;
import pty.smc.MapLeaves;

import ma.newick.NewickParser;
import nuts.math.Fct;
import nuts.math.Sampling;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.Tree;
import nuts.util.Arbre.ArbreMap;

/**
 * The size of the symmetric difference between 2 sets, i.e.
 * 
 * l(A,B) = |(A U B) \ (A V B)|
 * 
 * @author bouchard
 *
 * @param <T>
 */
public class SymmetricDiff<T> implements LossFct<Set<T>>
{
  public static <S> Set<S> filterLowCounts(Counter<S> counter, double threshold)
  {
    Set<S> result = new HashSet<S>();
    for (S key : counter.keySet())
      if (counter.getCount(key) > threshold)
        result.add(key);
    return result;
  }
  
  public static class ConsensusProcessor implements ParticleProcessor<PartialCoalescentState>
  {
    private Counter<Set<Taxon>> flatClades = new Counter<Set<Taxon>>();
    public void process(PartialCoalescentState state, double weight)
    {
      flatClades.incrementAll(state.allClades(), weight);
    }
    public Set<Set<Taxon>> consensus(double threshold)
    {
      if (threshold < 0.5) throw new RuntimeException();
      return filterLowCounts(flatClades, threshold);
    }
  }
  /**
   * @param <S>
   * @param unrootedTreeWithArbitraryRooting : should not be rooted at a leaf!!
   * @return
   */
  public static <S> Set<Set<S>> cladesFromUnrooted(Arbre<S> unrootedTreeWithArbitraryRooting)
  {
    Set<Set<S>> clades = clades(unrootedTreeWithArbitraryRooting);
    clades = addComplements(clades);
    clades.remove(new HashSet<Set<S>>());
    return clades;
  }
  /**
   * For each set, add the complement to the set of sets
   * Typical use: given an unrooted tree represented as a directed
   * tree with an arbitrary root (which is the way prescribed by newick),
   * complete the set of clades constraints extracted from clade(.)
   * @param original
   * @return
   */
  private static <S> Set<Set<S>> addComplements(Set<Set<S>> original)
  {
    Set<Set<S>> result = new HashSet<Set<S>>(original);
    Set<S> all = allLeaves(original);
    for (Set<S> clade : original)
    {
      HashSet<S> complement = new HashSet<S>(all);
      complement.removeAll(clade);
      result.add(complement);
    }
    return result;
  }
  
  public static <S> Set<S> allLeaves(Set<Set<S>> clades)
  {
    Set<S> result = new HashSet<S>();
    for (Set<S> clade : clades)
      result.addAll(clade);
    return result;
  }
  
  /**
   * Here, each elt in the symmetric difference is a clade,
   * where the clade is represented by a set of languages
   */
  public static final SymmetricDiff<Set<Taxon>> CLADE_SYMMETRIC_DIFFERENCE
    = new SymmetricDiff<Set<Taxon>>();
  
  public double loss(Set<T> s1, Set<T> s2)
  {
    return symmetricDifferenceSize(s1,s2);
  }
  public static <S> int symmetricDifferenceSize(Set<S> s1, Set<S> s2)
  {
    int result = 0;
    for (S elt : s1) 
      if (!s2.contains(elt)) 
        result++;
    for (S elt : s2) 
      if (!s1.contains(elt)) 
        result++;
//    System.out.println(result);
    return result;
  }
  
  public static int cladesSymmetricDifferenceSize(Arbre<Taxon> s1, Arbre<Taxon> s2, MapLeaves bijection)
  {
    return cladesSymmetricDifferenceSize(clades(s1),clades(s2),bijection);
  }
  public static int cladesSymmetricDifferenceSize(Set<Set<Taxon>> s1, Set<Set<Taxon>> s2, MapLeaves bijection)
  {
    Set<Set<Taxon>> mapped = bijection.mapClades(s1);
    Set<Set<Taxon>> filtered = bijection.filterClades(s2);
    return symmetricDifferenceSize(mapped,filtered);
  }
  
  public static ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> 
    createCladeProcessor()
  {
	  /*
	   * Inititalize the Decoder
	   * with a function<D,I> that converts the partial coalescent state D to 
	   * a set of clades I 
	   */
    return 
      new ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>>(
          new Fct<PartialCoalescentState, Set<Set<Taxon>>>() {
            public Set<Set<Taxon>> evalAt(PartialCoalescentState x)
            { return SymmetricDiff.clades(x.getUnlabeledArbre()); } });
  }
//  public static interface Bijection<S,T>
//  {
//    public Map<S,T> getForwardMap();
//    public Map<T,S> getBackwardMap();
//    public Set<S> getDomain();
//    public Set<T> getImage();
//    public Bijection<T,S> getInverse();
//  }
//  public static class HashBijection<S,T> implements Bijection<S,T>
//  {
//    private final Map<S,T> fwd;
//    private final Map<T,S> bwd;
//    public HashBijection(Map<S, T> fwd, Map<T, S> bwd)  
//    {
//      this.fwd = CollUtils.archive(fwd);
//      this.bwd = CollUtils.archive(bwd);
//    }
//    public Map<S,T> getForwardMap() { return fwd; }
//    public Map<T,S> getBackwardMap(){ return bwd; }
//    public Set<S> getDomain() { return getForwardMap().keySet(); }
//    public Set<T> getImage()  { return getBackwardMap().keySet();}
//    public Bijection<T,S> getInverse() { return new ReversedBijection(); }
//    private class ReversedBijection implements Bijection<T,S>
//    {
//      public Map<S,T> getBackwardMap() { return fwd; }
//      public Map<T, S> getForwardMap() { return bwd; }
//      public Set<T> getDomain() { return getForwardMap().keySet(); }
//      public Set<S> getImage()  { return getBackwardMap().keySet();}
//      public Bijection<S, T> getInverse() { return HashBijection.this; }
//    }
//  }
//  /**
//   * If x_0, x_1, .., x_n
//   * is a sequence of partial coalesc. states such that x_0 is initial and
//   * x_n is a full state, then:
//   * symmetricCladeDiff(clades(x_n), ref) = \sum_k^{n-1} ds(x_k, x_{k+1})
//   * where
//   * ds(x_k, x_{k+1}) = deltaSymmetricDiff(x_k, i, j, ref) : i,j is the pair coalesced in x_{k+1}
//   * @param state
//   * @param i
//   * @param j
//   * @param bijection a map from refs -> state
//   * @return
//   */
//  public static int deltaSymmetricDiff(PartialCoalescentState state, int i, int j, Set<Set<Language>> refs)
//  {
//    int result = 0;
//    if (state.nRoots() == 2) return 0; // root not counted
//    // 
//    Set<Language> newClade = state.mergedClade(i,j);
//    // check if the newly formed clade is in the reference
//    result += (refs.contains(newClade) ? 0 : 1);
//    // check if we violate for the first time some clade(s) in the reference
//    ref:for (Set<Language> ref : refs)
//      if (violates(ref, newClade))
//      {
//        // check not already accounted for
//        for (Set<Language> oldClade : state.allClades())
//          if (violates(ref, oldClade))
//            continue ref;
//        result++;
//      }
//    return result;
//  }  
//  /**
//   * If x_0, x_1, .., x_n
//   * is a sequence of partial coalesc. states such that x_0 is initial and
//   * x_n is a full state, then:
//   * symmetricCladeDiff(clades(x_n), ref) = \sum_k^{n-1} ds(x_k, x_{k+1})
//   * where
//   * ds(x_k, x_{k+1}) = deltaSymmetricDiff(x_k, i, j, ref) : i,j is the pair coalesced in x_{k+1}
//   * @param state
//   * @param i
//   * @param j
//   * @return
//   */
//  public static int deltaSymmetricDiff(PartialCoalescentState state, int i, int j, Set<Set<Language>> refs, MapLeaves ml)
//  {
//    int result = 0;
//    if (state.nRoots() == 2) return 0; // root not counted
//    // 
//    Set<Language> newClade = state.mergedClade(i,j);
//    // check if the newly formed clade is in the reference
//    result += (ml.contains(refs,newClade) ? 0 : 1);
//    // check if we violate for the first time some clade(s) in the reference
//    ref:for (Set<Language> ref : refs)
//      if (violates(ref, newClade, ml))
//      {
//        // check not already accounted for
//        for (Set<Language> oldClade : state.allClades())
//          if (violates(ref, oldClade, ml))
//            continue ref;
//        result++;
//      }
//    return result;
//  }
//  public static boolean violates(Set<Language> s1, Set<Language> s2, MapLeaves ml)
//  {
//    return ml.intersects(s1,s2) && !ml.containsAll(s2, s1) && !ml.containsAll(s1, s2);
//  }
  public static <T> boolean violates(Set<T> s1, Set<T> s2)
  {
    return CollUtils.intersects(s1,s2) && !s1.containsAll(s2) && !s2.containsAll(s1);
  }
  public static <T> boolean violatesOne(Set<T> s, Set<Set<T>> sets)
  {
    for (Set<T> s2 : sets)
      if (violates(s, s2))
        return true;
    return false;
  }
  public static int deltaSymmetricDiff(PartialCoalescentState state, int i, int j, 
      Set<Set<Taxon>> refs)
  {
    return deltaSymmetricDiff(state,i,j,refs,null);
  }
  public static int deltaSymmetricDiff(
      PartialCoalescentState state, 
      Set<Set<Taxon>> allCladesInState, 
      int i, int j, 
      Set<Set<Taxon>> preMappedRefs, 
      MapLeaves bijection)
  {
    int result = 0;
    if (state.nRoots() == 2) return 0; // root not counted
    Set<Taxon> newClade = state.mergedClade(i,j);
    if (bijection != null)
      newClade = bijection.restrict(newClade);
    if (!newClade.isEmpty()) // happens when there is a sparse bijection
    {
      // check if the newly formed clade is in the reference
      result += (preMappedRefs.contains(newClade) ? 0 : 1);
      // check if we violate for the first time some clade(s) in the reference
      ref:for (Set<Taxon> ref : preMappedRefs)
        if (violates(ref, newClade))
        {
          // check not already accounted for
          for (Set<Taxon> oldClade : allCladesInState)
            if (violates(ref, oldClade))
              continue ref;
          result++;
        }
    }
    return result;
  }
  public static int deltaSymmetricDiff(
      PartialCoalescentState state, 
      int i, int j, 
      Set<Set<Taxon>> refs, 
      MapLeaves bijection)
  {
    Set<Set<Taxon>> allCladesInState = state.allClades();
    if (bijection != null)
    {
      refs = bijection.mapClades(refs);
      allCladesInState = bijection.filterClades(allCladesInState);
    }
    return deltaSymmetricDiff(state, allCladesInState, i, j, refs, bijection);
//    
//    // 
//    if (!newClade.isEmpty()) // happens when there is a sparse bijection
//    {
//      Set<Set<Language>> allCladesInState = state.allClades();
//      if (bijection != null)
//        allCladesInState = bijection.filterClades(allCladesInState);
//      // check if the newly formed clade is in the reference
//      result += (refs.contains(newClade) ? 0 : 1);
//      // check if we violate for the first time some clade(s) in the reference
//      ref:for (Set<Language> ref : refs)
//        if (violates(ref, newClade))
//        {
//          // check not already accounted for
//          for (Set<Language> oldClade : allCladesInState)
//            if (violates(ref, oldClade))
//              continue ref;
//          result++;
//        }
//    }
//    return result;
  }
//  public static <S,T> Set<T> mapClade(Set<S> set, Map<S,T> map)
//  {
//    HashSet<T> result = new HashSet<T>();
//    for (S s : set)
//      if (map.keySet().contains(s))
//        result.add(map.get(s));
//    return result;
//  }

  
  public static void main(String [] args)
  {
    // load two, compute normalized
    NewickParser np;
    try {
      
      np = new NewickParser(IOUtils.openIn("data/generatedExperiment/tree.newick"));
      Tree<String> tree = np.parse();
      Arbre<String> ar = Arbre.tree2Arbre(tree);
      
      np = new NewickParser(IOUtils.openIn("data/generatedExperiment/contml50k.newick"));
      Tree<String> tree2 = np.parse();
      Arbre<String> ar2 = Arbre.tree2Arbre(tree2);

      System.out.println(normalizedSymmetricCladeDiff(ar,ar2));
    
    } catch (Exception e) { throw new RuntimeException(e); }
//    int nLeaves = 10;
//    Random rand = new Random();
//    for (int i = 0; i < 1000; i++)
//    {
//      PartialCoalescentState state = SymmetryTest.getInitState(nLeaves,5);
//      // generate a ref
//      PriorPriorKernel kernel = new PriorPriorKernel(state);
//      while (!state.isFinalState()) 
//        state = kernel.next(rand, state).getFirst();
//      Arbre<Language> _refs = state.getUnlabeledArbre();
//      Set<Set<Language>> refs = clades(_refs);
//      System.out.println("Ref:\n" + _refs.deepToString());
//      // generate another one and check deltas
//      int sum = 0;
//      state = SymmetryTest.getInitState(nLeaves,5);
//      while (!state.isFinalState())
//      {
//        // pick a pair to merge
//        List<Integer> sampledIndices = Sampling.sampleWithoutReplacement(rand, state.nRoots(), 2);
//        PartialCoalescentState next = state.coalesce(
//            sampledIndices.get(0), sampledIndices.get(1), 
//            rand.nextDouble());
//        int cur = deltaSymmetricDiff(state, sampledIndices.get(0), sampledIndices.get(1), refs);
//        System.out.print(cur + " ");
//        sum += cur;
//        state = next;
//      }
//      int dir = symmetricCladeDiff(_refs, state.getUnlabeledArbre());
//      System.out.println("\nOther:\n" + state.getUnlabeledArbre().deepToString());
//      System.out.println("Incremental:" + sum);
//      System.out.println("Direct:" + dir);
//      if (sum != dir) throw new RuntimeException();
//    }
  }
  
//  private static int _deltaSymmetricDiff(Set<Set<Language>> refs, Set<Language> leftClade, Set<Language> rightClade)
//  {
//    for (Set<Language> ref : refs) 
//      if (ref.containsAll(rightClade) && !CollUtils.intersects(ref,leftClade))
//        return 1;
//    return 0;
//  }
  
  /**
   * @return symmetricDifferenceSize(clades(t1),clades(t2))
   */
  public static <S> int symmetricCladeDiff(Arbre<S> t1, Arbre<S> t2)
  {
    // a set of partition (partitions are sets of sets)
    Set<Set<S>> clades1 = clades(t1),
                clades2 = clades(t2);
    return symmetricDifferenceSize(clades1, clades2);
  }
  public static <S> int maxSymmetricCladeDiff(Arbre<S> t1, Arbre<S> t2)
  {
    return t1.nodes().size() + t2.nodes().size();
  }
  public static <S> double normalizedSymmetricCladeDiff(Arbre<S> t1, Arbre<S> t2)
  {
    return ((double) symmetricCladeDiff(t1,t2))
                      /
           ((double) maxSymmetricCladeDiff(t1,t2));
  }
  /**
   * for all node in the rooted tree, its clade is the set of it desc. leaves 
   * @param <S>
   * @param t
   * @return the set of all clades
   */
  public static <S> Set<Set<S>> clades(Arbre<S> t)
  {
    Set<Set<S>> clades = new HashSet<Set<S>>();
    clades(leafSetArbre(t), clades);
    return clades;
  }
  private static <S> void clades(
      Arbre<Set<S>> t1, 
      Set<Set<S>> clades)
  {
    clades.add(t1.getContents());
    for (Arbre<Set<S>> child : t1.getChildren())
      clades(child, clades);
  }
  public static <S> Set<Set<S>> complete(Set<Set<S>> clades)
  {
    Set<Set<S>> result = new HashSet<Set<S>>(clades);
    Set<S> all = new HashSet<S>();
    for (Set<S> clade : clades)
      all.addAll(clade);
    result.add(all);
    for (S s : all)
      result.add(Collections.singleton(s));
    return result;
  }
  /**
   * reconstruct a tree given a set of clades
   * @param <S>
   * @param clades
   * @return
   */
  public static <S> Arbre<S> clades2arbre(Set<Set<S>> clades)
  {
    Arbre<S> result = _clades2arbre(complete(clades));
    return result;
  }
  private static <S> Arbre<S> _clades2arbre(Set<Set<S>> clades)
  {
    if (clades.size() == 1)
    {
      Set<S> singleton = clades.iterator().next();
      if (singleton.size() != 1) 
        throw new RuntimeException();
      return Arbre.arbre(singleton.iterator().next());
    }
    else
    {
      Set<Set<S>> remainder = removeLargest(clades);
      List<Arbre<S>> children = new ArrayList<Arbre<S>>();
      while (remainder.size() > 0)
      {
        Set<Set<S>> child = extractChild(remainder);
        children.add(_clades2arbre(child));
      }
      return Arbre.arbre(null,children);
    }
  }
  private static <S> Set<Set<S>> removeLargest(Set<Set<S>> rem)
  {
    Set<Set<S>> result = new HashSet<Set<S>>(rem);
    int max = Integer.MIN_VALUE;
    Set<S> argmax = null;
    for (Set<S> clade : rem)
      if (clade.size() > max)
      {
        max = clade.size();
        argmax = clade;
      }
    result.remove(argmax);
    return result;
  }
  private static <S> Set<Set<S>> extractChild(Set<Set<S>> rem)
  {
    // find the largest clade,
    int max = Integer.MIN_VALUE;
    Set<S> argmax = null;
    for (Set<S> clade : rem)
      if (clade.size() > max)
      {
        max = clade.size();
        argmax = clade;
      }
   // return it as well as all its child 
    Set<Set<S>> result = new HashSet<Set<S>>();
    // also remove them from rem
    Iterator<Set<S>> iter = rem.iterator();
    while (iter.hasNext())
    {
      Set<S> current = iter.next();
      if (argmax.containsAll(current))
      {
        result.add(current);
        iter.remove();
      }
    }
    return result;
  }
  
  /**
   * create a new Arbre where each node contains
   * the set of leaves under this leave
   * @param <T>
   * @param a
   * @return
   */
  public static <T> Arbre<Set<T>> leafSetArbre(Arbre<T> a)
  {
    return a.postOrderMap(new ArbreMap<T,Set<T>>() {
      @Override
      public Set<T> map(Arbre<T> currentDomainNode)
      {
        Set<T> result = new HashSet<T>(currentDomainNode.nLeaves());
        if (currentDomainNode.isLeaf())
          result.add(currentDomainNode.getContents());
        else
          for (Set<T> child : getChildImage())
            result.addAll(child);
        return result;
      }
    });
  }
}
