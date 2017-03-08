package pty;
import fig.basic.UnorderedPair;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import nuts.math.Sampling;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class RandomUnrootedTrees
{
  public static UnrootedTree sampleExponentialBranchesUniformTopology(Random rand, int nTaxa, double rate)
  {
    return sampleExponentialBranchesUniformTopology(rand, RandomRootedTrees.generateLeafNames(nTaxa), rate);
  }
  public static UnrootedTree sampleExponentialBranchesUniformTopology(Random rand, Collection<Taxon> taxa, double rate)
  {
    RootedTree coal = RandomRootedTrees.sampleCoalescent(rand, taxa, rate);
    UnrootedTree rooted = UnrootedTree.fromRooted(coal);
    Map<UnorderedPair<Taxon,Taxon>,Double> newBLs = map();
    for (UnorderedPair<Taxon,Taxon> key : rooted.edges())
      newBLs.put(key, Sampling.sampleExponential(rand, 1.0/rate));
    return new UnrootedTree(rooted.getTopology(), newBLs);
  }
}
