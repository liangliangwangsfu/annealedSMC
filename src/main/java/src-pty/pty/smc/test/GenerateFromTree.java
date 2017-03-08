package pty.smc.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.Taxon;
import ma.newick.NewickParser;
import ma.newick.ParseException;
import nuts.io.IO;
import nuts.util.Arbre;
import nuts.util.Tree;


public class GenerateFromTree implements Runnable {
    @Option public Random rand = new Random(1);
    @Option public double scale = 0.001;
    @Option public double root = 0.5;
    @Option public int nsites = 10000;
    @Option(required=true) public String treeFile = "";
    public int fixed = 0;

	public File generate () {
	  LogInfo.logsForce("Scale=" + scale);
	  LogInfo.logsForce("Root="+root);
	  LogInfo.logsForce("nSites="+nsites);
	  LogInfo.logsForce("Tree=" + treeFile);
		NewickParser np;
		try {
			np = new NewickParser(IOUtils.openIn(treeFile));
			Tree<String> tree = np.parse();
			Arbre<String> ar = Arbre.tree2Arbre(tree);
			HashMap<Arbre<String>,double[]>  node2sequenceMap =
				 new HashMap<Arbre<String>, double[]>();
			Map<Taxon,Double> branchLengths = np.getBranchLengths();
			LogInfo.logsForce(ar.deepToString());
			
//		    BufferedWriter bw = new BufferedWriter (new FileWriter("simulated.txt"));
			File result = new File(Execution.getFile("simulated.txt"));
			PrintWriter bw = IOUtils.openOutHard(result);
			
			List<Arbre<String>> preorder = ar.nodes();
			for (Arbre<String> node:preorder){
				LogInfo.logsForce ("Examining node:" + node.getContents());
				double[] tmp = new double[nsites];
				if (node.isRoot()) {
				    for (int i = 0; i < nsites; i++)
				    	tmp[i] = root;  
				} else {
					double[] parentsequence = node2sequenceMap.get(node.getParent());

					double bl = branchLengths.get(new Taxon(node.getContents()));
					LogInfo.logsForce ("Branch length = " + bl);
					LogInfo.logsForce("Parent = " + node.getParent());
					for (int i = 0; i < nsites; i++)
						tmp[i] = fix(parentsequence[i] + rand.nextGaussian()* Math.sqrt(scale * bl));
					
				}
				
				if (node.isLeaf()) {
					bw.write (node.getContents()+"\t");
					for (int i = 0; i < nsites; i++)
						bw.write (tmp[i]+"\t");
					bw.write ("\n");
				}
				node2sequenceMap.put(node, tmp );
			}
			bw.close ();
			LogInfo.logsForce ("Fixed " + fixed + " sites");
			return result;
			
		} catch (Exception e) { throw new RuntimeException(e); }
		
	}
	
	  private double fix(double d)
	  {
	    if (d < 0)
	    {
	    //  System.err.println("Less than zero");
	    	fixed++;
	      return 0;
	    }
	    else if (d > 1)
	    {
	    //  System.err.println("More than one");
	    	fixed++;
	      return 1;
	    }
	    else return d;
	  }

	  
	  public static void printArray (double[] x) {
		  for (int i = 0; i < x.length; i++)
			  System.out.print(x[i] +" \t");
		  System.out.print("\n");
	  }
	  
	  
	  
	  public static void main (String args[]) {
	    IO.run(args, new GenerateFromTree());
	  }

    @Override
    public void run()
    {
      generate();
    }
}
