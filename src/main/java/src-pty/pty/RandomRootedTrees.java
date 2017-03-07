package pty;
import fig.prob.SampleUtils;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import nuts.math.Sampling;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import nuts.util.MathUtils;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class RandomRootedTrees
{
  public static RootedTree sampleUnif(Random rand, int nTax)
  {
    return sampleUnif(rand, generateLeafNames(nTax));
  }
  
  public static RootedTree sampleUnif(Random rand, Collection<Taxon> leaves)
  { 
    return sample(rand, leaves, Double.NaN, false, true);
  }
  
  
  
  public static RootedTree sampleCoalescent(Random rand, int nTax, double rate)
  {
    return sampleCoalescent(rand, generateLeafNames(nTax), rate);
  }
  
  public static RootedTree sampleCoalescent(Random rand, Collection<Taxon> leaves, double rate)
  { 
    return sample(rand, leaves, rate, false);
  }
  
  public static RootedTree sampleYule(Random rand, int nTax, double rate)
  {
    return sampleYule(rand, generateLeafNames(nTax), rate);
  }
  
  public static RootedTree sampleYule(Random rand, Collection<Taxon> leaves, double rate)
  { 
    return sample(rand, leaves, rate, true);
  }
  private static RootedTree sample(Random rand, Collection<Taxon> leaves, double rate, boolean yule)
  { 
    return sample(rand, leaves, rate, yule, false);
  }
  private static RootedTree sample(Random rand, Collection<Taxon> leaves, double rate, boolean yule, boolean unif)
  { 
    if (unif && yule) throw new RuntimeException();
    List<RootedTree> current = list();
    List<Double> heights = list();
    double currentMaxHeight = 0.0;
    for (Taxon t : leaves)
    {
      current.add(RootedTree.Util.singleton(t));
      heights.add(0.0);
    }
    int i = 0;
    while (current.size() > 1)
    {
      // pick a delta
      final int n = current.size();
      final double currentRate = rate * (yule ? n : MathUtils.nChoose2(n));
      final double delta = unif ? Sampling.nextDouble(rand, 0.0, 1.0) : Sampling.sampleExponential(rand, 1.0/currentRate);
      // pop a pair at random to merge:
      List<Integer> pair = Sampling.sampleWithoutReplacement(rand, current.size(),2);
      RootedTree first = current.get(pair.get(0)), second = current.get(pair.get(1));
      double firstH = heights.get(pair.get(0)), secondH = heights.get(pair.get(1));
      current.remove(Math.max(pair.get(0),pair.get(1)));
      current.remove(Math.min(pair.get(0),pair.get(1)));
      heights.remove(Math.max(pair.get(0),pair.get(1)));
      heights.remove(Math.min(pair.get(0),pair.get(1)));
      final double newMaxHeight = currentMaxHeight + delta;
      currentMaxHeight = newMaxHeight;
      double left = currentMaxHeight - firstH;
      double right = currentMaxHeight - secondH;
      RootedTree newTree = RootedTree.Util.coalesce(new Taxon("internal_" + (i++)), first, second, left, right);
      current.add(newTree);
      heights.add(currentMaxHeight);
    }
    return current.get(0);
  }
  
  public static List<Taxon> generateLeafNames(int nTax)
  {
    List<Taxon> rst = list();
    for (int i =0 ; i < nTax ;i++)
      rst.add(new Taxon("leaf_" + i));
    return rst;
  }

  public static void main(String [] args)
  {
    Random rand = new Random(1);
    RootedTree rt = (sampleYule(rand, 10, 1.0));
    System.out.println(rt);
  }
}
