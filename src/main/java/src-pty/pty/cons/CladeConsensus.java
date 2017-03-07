package pty.cons;
import java.io.*;
import java.util.*;

import fig.basic.UnorderedPair;
import goblin.Taxon;

import pty.RootedTree;
import pty.UnrootedTree;
import pty.eval.SymmetricDiff;
import nuts.util.Arbre.ArbreMap;
import nuts.util.CollUtils.*;
import nuts.util.Arbre;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class CladeConsensus
{
  public static Arbre<Taxon> consensusTopology(Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> counter, double threshold)
  {
    if (threshold < 0.5)
      throw new RuntimeException();
    Set<Set<Taxon>> rootedClades = set();
    UnorderedPair<Set<Taxon>,Set<Taxon>> aBipart = pick(counter.keySet());
    Taxon aTaxon = pick(aBipart.getFirst());
    for (UnorderedPair<Set<Taxon>, Set<Taxon>> key : counter.keySet())
      if (counter.getCount(key) >= threshold)
        rootedClades.add(rootClade(key, aTaxon));
    
    for (Taxon t : aBipart.getFirst()) rootedClades.add(Collections.singleton(t));
    for (Taxon t : aBipart.getSecond()) rootedClades.add(Collections.singleton(t));
    rootedClades.add(union(aBipart.getFirst(),aBipart.getSecond()));

    Arbre<Taxon> a =  SymmetricDiff.clades2arbre(rootedClades);
    return a.preOrderMap(new ArbreMap<Taxon,Taxon>() {
      private int n = 0;
      @Override
      public Taxon map(Arbre<Taxon> currentDomainNode)
      {
        if (currentDomainNode.getContents() == null)
          return new Taxon("internal_" + (n++));
        else
          return currentDomainNode.getContents();
      }
    });
  }
  
  public static UnrootedTree consensusUnrootedTree(Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> counter, double threshold)
  {
    if (threshold < 0.5 || threshold > 1) throw new RuntimeException();
    Arbre<Taxon> topology = consensusTopology(counter, threshold);
    Map<Arbre<Taxon>,Set<Taxon>> leaves = Arbre.leavesMap(topology);
    Set<Taxon> allTaxa = leaves.get(topology);
    Map<Taxon,Double> bls = map();
    for (Arbre<Taxon> subt : topology.nodes())
      if (!subt.isRoot())
      {
        // create the key:
        Set<Taxon> set1 = leaves.get(subt);
        Set<Taxon> complement = set(allTaxa);
        complement.removeAll(set1);
        UnorderedPair<Set<Taxon>, Set<Taxon>> key = new UnorderedPair<Set<Taxon>, Set<Taxon>>(set1, complement);
        double currentBL = counter.getCount(key);
        bls.put(subt.getContents(), currentBL);
      }
    RootedTree rt = new RootedTree.Util.RootedTreeImpl(topology, bls);
    return UnrootedTree.fromRooted(rt);
  }
  
  public static void main(String [] args)
  {
    // load 
//    UnrootedTree temp = UnrootedTree.fromNewick(new File("scratch/toy.txt"));
    UnrootedTree temp = UnrootedTree.fromNewick(new File("e/575.exec/WALS-chain-0/mle.newick"));
    System.out.println(temp);
    Counter<UnorderedPair<Set<Taxon>, Set<Taxon>>> counter = temp.inducedBiPartitions2BranchMap();
    for (UnorderedPair<Set<Taxon>, Set<Taxon>> key : set(counter.keySet()))
      counter.setCount(key, 1.0);
    UnrootedTree rt = consensusUnrootedTree(counter, 0.5);
    System.out.println(rt);
  }

  private static Set<Taxon> rootClade(UnorderedPair<Set<Taxon>, Set<Taxon>> bipartition, Taxon out)
  {
    if (bipartition.getFirst().contains(out)) 
      return bipartition.getSecond();
    else if (bipartition.getSecond().contains(out))
      return bipartition.getFirst();
    else
      throw new RuntimeException();
  }
}
