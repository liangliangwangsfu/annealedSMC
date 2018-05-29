package smcsampler;

import static nuts.io.IO.i;

import java.io.File;
import java.io.PrintWriter;
import ev.poi.processors.TreeDistancesProcessor;
import ev.to.NJ;
import fig.basic.IOUtils;
import fig.basic.Option;
import fig.exec.Execution;
import nuts.io.CSV;
import nuts.io.IO;
import pty.RootedTree;
import pty.UnrootedTree;
import pty.RootedTree.RootedTreeProcessor;
import pty.io.TreeEvaluator;
import pty.io.TreeEvaluator.TreeMetric;

public class ConcensusTreeFromTreeList  implements Runnable{
	@Option public String treeListFile=null;
	@Option public String outConcensusTreeFile=null;
	@Option public String refTree=null;
	@Option public String outTreeDistanceFile=null;
	
	
	public UnrootedTree getConcensusTree()
	{		
		TreeDistancesProcessor tdp = new TreeDistancesProcessor();		
		 processTreeList(new File(treeListFile), tdp);
		//MrBayes.processMrBayesTrees(new File(treeListFile), tdp);
		return tdp.getConsensus();		
	}
	
	
	public static void processTreeList(File file, RootedTreeProcessor rtp)
	{
		for (String line : i(file))
		{
			if(line.length()>0)
				rtp.process(RootedTree.Util.fromNewickString(line));
		}
	}
	
	public static void computeTreeMetrics(UnrootedTree inferred, String refTree, String treeListFile, String filename)
	{				
		UnrootedTree goldut = UnrootedTree.fromRooted(RootedTree.Util.fromNewickString(IO.f2s(refTree)));
		PrintWriter out = IOUtils.openOutEasy(new File(Execution.getFile(filename)));
		out.println(CSV.header("Metric", "Value", "treeList", "refTree"));
		for (TreeMetric tm : TreeEvaluator.coreTreeMetrics)
		{
			double value = tm.score(inferred, goldut);
			out.println(CSV.body(tm.toString(), value,  treeListFile, refTree));
		}
		out.close();		
	}
	
	
	@Override
	public void run() {
		UnrootedTree inferred = getConcensusTree();		
		IO.writeToDisk(new File(Execution.getFile(outConcensusTreeFile)), getConcensusTree().toNewick());
		computeTreeMetrics(inferred, refTree,  treeListFile, outTreeDistanceFile);		
		
	}
	
	public static void main(String[] args)
	{
		IO.run(args, new ConcensusTreeFromTreeList(),"nj", NJ.class);
	}

}
