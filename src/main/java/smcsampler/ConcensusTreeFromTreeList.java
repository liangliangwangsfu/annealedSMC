package smcsampler;

import java.io.File;
import ev.poi.processors.TreeDistancesProcessor;
import ev.to.NJ;
import fig.basic.Option;
import nuts.io.IO;
import pty.UnrootedTree;


public class ConcensusTreeFromTreeList  implements Runnable{
	@Option public String treeListFile=null;
	@Option public String outConcensusTreeFile=null;
	
	public UnrootedTree getConcensusTree()
	{		
		TreeDistancesProcessor tdp = new TreeDistancesProcessor();
		MrBayes.processMrBayesTrees(new File(treeListFile), tdp);
		return tdp.getConsensus();		
	}
	
	@Override
	public void run() {
		IO.writeToDisk(new File(outConcensusTreeFile), getConcensusTree().toNewick());
	}
	
	public static void main(String[] args)
	{
		IO.run(args, new ConcensusTreeFromTreeList(),"nj", NJ.class);
	}

}
