package pty.eval;
import java.io.*;
import java.util.*;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import pty.UnrootedTree;

import nuts.io.IO;
import nuts.lang.CollectionUtils;
import nuts.math.EqClasses;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Tree;

import fig.basic.IOUtils;
import fig.basic.UnorderedPair;
import goblin.Taxon;

import ma.newick.NewickParser;
import ma.newick.ParseException;

public class Purity
{
  /**
   * See footnote 2 of the Heller et al., Bayesian Hierarchical Clustering
   * 
   * @param <S> type of the species
   * @param <T> type of the species clusters
   * @param tree
   * @param labels a map from species to their cluster
   * @return
   */
  public static <S,T> double purity(Tree<S> tree, Map<S,T> labels)
  {
    return purity(SymmetricDiff.clades(Arbre.tree2Arbre(tree)), labels);
  }
//  {
//    Map<UnorderedPair<S,S>,Tree<S>> lcas = lcas(tree);
////    if (labels.size() > tree.getYield().size())
////      System.err.println("Warning: more labels than leaves in WalsEval.purity()");
//    // take a leave u.a.r.
//    SummaryStatistics mainAvg = new SummaryStatistics();
//    for (S leave : tree.getYield())
//    {
//      // take another leave u.a.r. among those with same label
//      SummaryStatistics subAvg = new SummaryStatistics();
//      T cLabel = labels.get(leave);
//      for (S other : tree.getYield())
//        if (!other.equals(leave) && cLabel.equals(labels.get(other)))
//          subAvg.addValue(
//             ratio(
//                 lcas.get(new UnorderedPair<S,S>(leave,other)), 
//                 labels, 
//                 cLabel));
//      if (subAvg.getN() > 0)
//        mainAvg.addValue(subAvg.getMean());
//    }
//    return mainAvg.getMean();
//  }
  public static <T> double purity(UnrootedTree t, Map<Taxon, T> labels)
  {
    return purity(t.clades(), labels);
  }
  public static <T> double recall(Set<T> reference, Set<T> guess)
  {
    return CollUtils.inter(reference, guess).size() / reference.size();
  }
  public static <S,T> double purity(Set<Set<S>> clades, Map<S,T> labels)
  {
    Set<S> leaves = leaves(clades);
//    if (labels.size() > tree.getYield().size())
//      System.err.println("Warning: more labels than leaves in WalsEval.purity()");
    // take a leave u.a.r.
    SummaryStatistics mainAvg = new SummaryStatistics();
    for (S leave : leaves)
    {
      // take another leave u.a.r. among those with same label
      SummaryStatistics subAvg = new SummaryStatistics();
      T cLabel = labels.get(leave);
      if (cLabel != null)
        for (S other : leaves)
        {
          T otherLabel =  labels.get(other);
          if (!other.equals(leave) && otherLabel != null && cLabel.equals(otherLabel)) //labels.get(other)))
            subAvg.addValue(
               ratio(
                   lca(leave,other,clades), 
                   labels, 
                   cLabel));
        }
      if (subAvg.getN() > 0)
        mainAvg.addValue(subAvg.getMean());
    }
    return mainAvg.getMean();
  }
//  private static <S,T> Object safeGet(Map<S,T> labels, S key)
//  {
//    Object result = labels.get(key);
//    if (result == null)
//      result = key;
//    return result;
//  }
  /**
   * The smallest clade that contains both leaves
   */
  private static <S> Set<S> lca(final S leave1, final S leave2, final Set<Set<S>> clades)
  {
    Set<S> result = null;
    for (Set<S> clade : clades)
      if ((result == null || clade.size() < result.size()) &&
          clade.contains(leave1) && clade.contains(leave2))
        result = clade;
    return result;
  }
  private static <S> Set<S> leaves(Set<Set<S>> clades)
  {
    Set<S> result = new HashSet<S>();
    for (Set<S> clade : clades)
      if (clade.size() == 1)
        result.addAll(clade);
    return result;
  }
  
  public static <S,T> Map<T,Set<S>> partitionsUsedForEval(Tree<S> t, Map<S,T> labels)
  {
    Set<S> leaves = new HashSet<S>(t.getYield());
    Map<S,T> labelsCopy = new HashMap<S,T>(labels);
    labelsCopy.keySet().retainAll(leaves);
    Map<T,Set<S>> inverse = CollUtils.invert(labelsCopy);
    Iterator<T> iter = inverse.keySet().iterator();
    while (iter.hasNext())
      if (inverse.get(iter.next()).size() == 1)
        iter.remove();
    return inverse; 
  }
  
  public static void main(String [] args) throws IOException, ParseException
  {
    if (args.length != 2)
    {
      System.err.println("Evaluates the purity of a newick linguistic tree against");
      System.err.println("genus annotation in the file language.tab in wals_data");
      System.err.println("<newick-format-tree> <path-to-wals-languages-file>");
      return;
    }
    
    // read tree file
    NewickParser np = new NewickParser(IOUtils.openIn(args[0]));
    Tree<String> tree = np.parse();
    Map<String,String> lang2genus = new HashMap<String,String>();
    
    // read genus/family file
    for (String line : IO.i(args[1]))
      if (!line.startsWith("wals"))
      {
        String [] fields = line.split("\\t");
        lang2genus.put(fields[0], fields[4]);
      }
    
    // score
    System.out.println("Purity: " + purity(tree,lang2genus));
  }
  


  
  /**
   * @param <S>
   * @param <T>
   * @param subtree
   * @param labels
   * @return the ratio of leaves in a subtree that have the given label
   */
  private static <S,T> double ratio(Set<S> clade, Map<S,T> labels, T label)
  {
    SummaryStatistics ratio = new SummaryStatistics();
    for (S leaf : clade)
    {
      T leafLabel = labels.get(leaf);
      if (leafLabel == null) 
        ratio.addValue(0.0);
      else 
        ratio.addValue(leafLabel.equals(label) ? 1.0 : 0.0);
    }
    if (Double.isNaN(ratio.getMean()))
      throw new RuntimeException();
    return ratio.getMean();
  }
  
  public static <S> Map<UnorderedPair<S,S>,Tree<S>> lcas(Tree<S> tree)
  {
    Map<UnorderedPair<S,S>,Tree<S>> result = new HashMap<UnorderedPair<S,S>,Tree<S>>();
    _lcas(tree, result);
    return result;
  }
  private static <S> void _lcas(Tree<S> tree, Map<UnorderedPair<S,S>,Tree<S>> result)
  {
    for (int i = 0; i < tree.getChildren().size(); i++)
      for (int j = i+1; j < tree.getChildren().size(); j++)
        for (S first : tree.getChildren().get(i).getYield())
          for (S second : tree.getChildren().get(j).getYield())
            result.put(new UnorderedPair<S,S>(first,second),tree);
    for (Tree<S> child : tree.getChildren())
      _lcas(child,result);
  }
}
