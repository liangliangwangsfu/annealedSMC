package pty.io;
import java.io.File;
import java.util.Arrays;
import java.util.List;

import nuts.io.IO;
import nuts.math.Evaluator;
import nuts.math.Evaluator.EvaluationScorer;
import nuts.math.Evaluator.Parser;
import pty.RootedTree;
import pty.UnrootedTree;
import fig.basic.Option;


public class TreeEvaluator
{
  public static void main(String [] args)
  {
    IO.run(args, new TreeEvaluatorMain());
  }
  
  public static class TreeEvaluatorMain implements Runnable
  {
    @Option(required=true) public File refDirectory;
    @Option(required=true) public File guessDirectory;
    @Option public String refSuffix = null;
    @Option public String guessSuffix = null;
    @Override
    public void run()
    {
      Evaluator.evaluate(refDirectory, refSuffix, guessDirectory, guessSuffix, new RootedTreeParser(), treeEvaluationMetrics);
    }
  }
  
  public static class RootedTreeParser implements Parser<RootedTree>
  {
    @Override public RootedTree parse(File f) { return RootedTree.Util.load(f); }
  }
  
  public static abstract class TreeMetric implements EvaluationScorer<RootedTree>
  {
    public double worstScore() { return Double.POSITIVE_INFINITY; }
    public double score(UnrootedTree ut1, UnrootedTree ut2) { throw new RuntimeException(); }
    @Override public double score(RootedTree rt1, RootedTree rt2)
    {
      return score(UnrootedTree.fromRooted(rt1), UnrootedTree.fromRooted(rt2));
    }

		public double score(List<RootedTree> rt1, List<RootedTree> rt2) {
			double sum=0;
			for(int i=0;i<rt1.size();i++)
			sum+=score(UnrootedTree.fromRooted(rt1.get(i)),
					UnrootedTree.fromRooted(rt2.get(i)));			
			return sum;
		}

    @Override public boolean equals(Object obj) { return toString().equals(obj); }
    @Override public int hashCode() { return toString().hashCode(); }
  }
  @SuppressWarnings("unchecked")
  public static List<EvaluationScorer<RootedTree>> treeEvaluationMetrics = (List)Arrays.asList(
      new PartitionMetric(),
      new NormalizedPartitionMetric(),
      new TightlyNormalizedPartitionMetric(),
      new RobinsonFouldsMetric(),
      new KuhnerFelsenstein(),
      new NormalizedRobinsonFouldsMetric(),
      new NormalizedKuhnerFelsenstein()
      );
  
  public static List<TreeMetric> coreTreeMetrics = (List)Arrays.asList(
      new PartitionMetric(),
      new TightlyNormalizedPartitionMetric(),
      new RobinsonFouldsMetric(),
      new NormalizedRobinsonFouldsMetric(),
      new KuhnerFelsenstein()
      );
  
  public static class PartitionMetric extends TreeMetric
  {
    @Override public double score(UnrootedTree ref, UnrootedTree guess)  {  return UnrootedTree.partitionMetric(ref, guess); }
    @Override public String toString() { return "PartitionMetric"; }
  }
  public static class NormalizedPartitionMetric extends TreeMetric
  {
    @Override public double score(UnrootedTree ref, UnrootedTree guess)  {  return UnrootedTree.normalizedPartitionMetric(ref, guess); }
    @Override public String toString() { return "NormalizedPartitionMetric"; }
  }
  public static class TightlyNormalizedPartitionMetric extends TreeMetric
  {
    @Override public double score(UnrootedTree ref, UnrootedTree guess)  {  return UnrootedTree.tightlyNormalizedPartitionMetric(ref, guess); }
    @Override public String toString() { return "TightlyNormalizedPartitionMetric"; }
  }
  public static class RobinsonFouldsMetric extends TreeMetric
  {
    @Override public double score(UnrootedTree ref, UnrootedTree guess)  {  return UnrootedTree.robinsonFouldsMetric(ref, guess); }
    @Override public String toString() { return "RobinsonFouldsMetric"; }
  }

	public static class MatchingMetric extends TreeMetric {
		@Override
		public double score(UnrootedTree ref, UnrootedTree guess) {
			return UnrootedTree.matchingMetric(ref, guess);
		}

		@Override
		public String toString() {
			return "MatchingMetric";
		}
	}
  public static class KuhnerFelsenstein extends TreeMetric
  {
    @Override public double score(UnrootedTree ref, UnrootedTree guess)  {  return UnrootedTree.kuhnerFelsenstein(ref, guess); }
    @Override public String toString() { return "KuhnerFelsenstein"; }
  }
  public static class NormalizedRobinsonFouldsMetric extends TreeMetric
  {
    @Override public double score(UnrootedTree ref, UnrootedTree guess)  
    {  
      return UnrootedTree.normalizedRobinsonFouldsMetric(ref, guess); 
    }
    @Override public String toString() { return "NormalizedRobinsonFouldsMetric"; }
  }
  public static class NormalizedKuhnerFelsenstein extends TreeMetric
  {
    @Override public double score(UnrootedTree ref, UnrootedTree guess)  
    {  
      return UnrootedTree.normalizedKuhnerFelsenstein(ref, guess); 
    }
    @Override public String toString() { return "NormalizedKuhnerFelsenstein"; }
  }
}
