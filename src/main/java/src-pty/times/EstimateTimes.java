package times;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

//import algorithms.Barrier;
//import algorithms.MatrixPrint;

import Jama.Matrix;

import nuts.util.Arbre;
import nuts.util.Tree;

import fig.basic.IOUtils;
import goblin.Taxon;


import ma.newick.NewickParser;

public class EstimateTimes {
//	HashMap <String, Integer> popSizeMap;
//	Matrix a;
//	Matrix b;
//	Matrix c;
//	
//	public HashMap<String,Integer> readMap (String file) {
//		HashMap <String, Integer> tmp = new HashMap <String, Integer> ();
//		
//		try {
//			BufferedReader br = new BufferedReader (new FileReader (file));
//			String s;
//			while  ( (s = br.readLine()) !=null) {
//				StringTokenizer tok = new StringTokenizer (s);
//				String name = tok.nextToken();
//				int size = Integer.parseInt(tok.nextToken());
//				tmp.put(new String(name), size);
//				
//			}
//		}catch (Exception e) {
//			e.printStackTrace();
//		}
//		return tmp;
//	}
//	
//	public EstimateTimes (String file, String popSizeFile){
//		double[][] a;
//		double[] b;
//		double[][] c;
//		double[] init;
//		
//		popSizeMap = readMap (popSizeFile);
//		int numKnownPopSizes = 0 ;
//		
//		try {
//			NewickParser np = new NewickParser(IOUtils.openIn (file) );
//			Tree<String> tree= np.parse();
//			Arbre<String> ar = Arbre.tree2Arbre(tree);
//			
//			BufferedWriter bw = new BufferedWriter (new FileWriter ("tree.txt"));
//			bw.write (ar.deepToString());
//			bw.close ();
//			
//			Map<Taxon,Double> tmpbranchLengths = np.getBranchLengths();
//			Map<String, Double> branchLengths = new HashMap<String,Double> ();
//			for (Taxon l:tmpbranchLengths.keySet()) {
//				branchLengths.put(l.toString(), tmpbranchLengths.get(l));
//			}
//			
//			List<Arbre<String>> preorder = ar.nodes();
//			List<Arbre<String>> postorder = ar.nodesInPostOrder();
//			
//			HashMap <String, Integer> tmap = new HashMap<String, Integer> ();
//			HashMap <String, Integer> nmap = new HashMap <String, Integer> ();
//			
//			int numLeaves = ar.nLeaves();
//			int numInodes = preorder.size () - numLeaves;
//			int tmapcounter = 0;
//			int nmapcounter = numInodes;
//			
//			BufferedWriter bw1 =  new BufferedWriter (new FileWriter ("vars.txt"));
//			String var1 = "";
//			String var2 =  "";
//			int j  =0 ;
//			for (Arbre<String> node:preorder){
//				if (!node.isLeaf() ) { 
//					tmap.put(node.getContents(), tmapcounter);
//					var1 += tmapcounter + "\tTime_"+node.getContents()+"\n";
//					tmapcounter++;
//
//					if (!node.isRoot()) {
//						nmap.put(node.getContents(), nmapcounter);
//						var2 += nmapcounter + "\tN_"+node.getContents()+"\n";
//						nmapcounter ++;
//					}
//				} else {
//					j ++;
//					if ( !popSizeMap.containsKey(new String(node.getContents()))) {
//						nmap.put(node.getContents(), nmapcounter);
//						var2 += nmapcounter + "\tN_"+node.getContents()+"\n";
//						nmapcounter ++;
//					} else 
//						numKnownPopSizes ++;
//				}								
//			}
//			
//			System.out.println ("sizse = " + preorder.size() + "\t" +  j);
//			System.out.println ("numLeaves = "  + numLeaves);
//			System.out.println ("tmapcounter = " + tmapcounter);
//			System.out.println ("nmapcounter = " + nmapcounter);
//			
//			bw1.write(var1);
//			bw1.write(var2);
//			bw1.close ();
//			
//			int nvar = nmapcounter;
//			int nconstraints = numLeaves + numInodes  - 1;
//			a = new double[nconstraints][nvar];
//			b = new double[nconstraints];
//			c = new double[nvar + numInodes - 1][nvar];
//			init = new double[nvar];
//			
//			int acounter =  0;
//			int ccounter = 0;
//			
//			for ( ccounter = 0 ; ccounter < nvar; ccounter++)
//				c[ccounter][ccounter] = -1;
//			
//
//			
//			for (Arbre<String> node:preorder){
//				if (node.isLeaf() ) {
//					Arbre<String> par = node.getParent();
//					int t1 = tmap.get(par.getContents());
//					double bl = branchLengths.get (node.getContents());
//					
//					if ( popSizeMap.containsKey(node.getContents())) {
//						int n1 = popSizeMap.get(node.getContents());					
//						a[acounter][t1] = 1;
//						b[acounter] = n1*bl;
//						acounter ++;
//					} else {
//						int n = nmap.get (node.getContents());
//						a[acounter][t1] = 1;
//						a[acounter][n] = -bl;
//						acounter ++;
//					}
//				} else if (node.isRoot()) {					
//				} else {
//					Arbre<String> par = node.getParent();
//					int t2 = tmap.get(node.getContents());
//					int t1 = tmap.get(par.getContents());
//					int n1 = nmap.get(node.getContents());
//					double bl = branchLengths.get (node.getContents());
//					a[acounter][t1] = 1;
//					a[acounter][t2] = -1;
//					a[acounter][n1] = -bl;
//					b[acounter] = 0;
//					acounter ++;
//					
//					c[ccounter][t1] = -1;
//					c[ccounter][t2] = 1;
//					ccounter++;					
//				}									
//			}
//			
//			this.a = new Matrix (a,nconstraints,nvar);
//			this.b = new Matrix (b,nconstraints);
//			this.c = new Matrix (c, nvar + numInodes - 1, nvar);
//			
//			for (Arbre<String> node:postorder){
//				if (node.isLeaf()) {
//					
//				} else {
//					
//				}
//			}			
//			
//		} catch (Exception e){
//			e.printStackTrace();
//		}
//	}
//	
//	public void printMatrices  () {
//		MatrixPrint.printMatln("a =",a);
//		MatrixPrint.printMatln("b =",b);
//		MatrixPrint.printMatln("c =",c);
//	}
//	
//	public void printMFile () {
//		MatrixPrint.printFile("a.txt", a);
//		MatrixPrint.printFile("b.txt", b);
//		MatrixPrint.printFile("c.txt", c);
//		
//		String file = "";
//		file = "a = load (\'a.txt\');\n";
//		file += "b = load (\'b.txt\');\n";
//		file += "c = load (\'c.txt\');\n";
//		file += "n = size(a,2);\n";
//		file += "cvx_begin;\n";
//		file += "variable x(n);\n";
//		file += "minimize norm (a*x-b);\n";
//		file += "subject to \n";
//		file += "c * x <= 0; \n";		
//		file += "cvx_end\n";
//		file += "x = x*20;\n";
//		file += "save -ASCII \'x.txt\' \'x\';\n";
//		
//		try {
//			FileWriter fw = new FileWriter ("estimatetimes.m");
//			fw.write(file);
//			fw.close ();
//			
//			Runtime run = Runtime.getRuntime() ;
//			Process pr = run.exec("./runmatlab.sh") ;
//			pr.waitFor();
//			
//		} catch (Exception e){
//			e.printStackTrace();
//		}
//	}
//
////	
////	public Matrix solve () { 
////		Barrier b  = new Barrier ()
////	}
////	
//	public static void main (String args[]) {
//		EstimateTimes et = new EstimateTimes (args[0], args[1]);
//		et.printMFile();
//	}
}