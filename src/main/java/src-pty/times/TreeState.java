package times;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.StringTokenizer;

import pty.UnrootedTree;

import Jama.Matrix;

import ma.newick.NewickParser;
import nuts.io.IO;
import nuts.util.Arbre;
import nuts.util.Tree;
import fig.basic.IOUtils;
import fig.basic.Pair;


import goblin.Taxon;

public class TreeState implements State{
	
	
	private HashMap<Set<String>, String > cladesMap;
	private HashMap<String, Arbre<String> > arbreMap;
	private HashMap<String, Double> branchLengths;
	private Arbre<String> ar;
	private TimeState timeState;
	private RateState rateState;
	public static double alpha = 0.05;
	
	public TimeState getTimeState () {
		return timeState;
	}
	
	public RateState getRateState () {
		return rateState;
	}
	
	public double getLogDensity() {
		// TODO Auto-generated method stub
		double logd = jacobian ();
		logd += timeState.getLogDensity() + rateState.getLogDensity();
		return logd;
	}	

	
	private double jacobian () {
		int n = ar.nLeaves();
		n = (n - 1)  + ( 2 * n - 2);
		Matrix m = new Matrix (n,n);
		
		HashMap <String, Integer> time2index =  new HashMap <String, Integer> ();
		HashMap <String, Integer> edge2index = new  HashMap <String, Integer> ();
		int index =  0 ;
		for (Arbre<String> a:ar.nodes()) {		
			if ( !a.isLeaf() ) 
				time2index.put (a.getContents(),index++);
		}
		
		for (String s:branchLengths.keySet()) {
			edge2index.put(s, index++);
		}
		
		for (String s:branchLengths.keySet()) {
			double r = rateState.getRate(s);
			double v = branchLengths.get(s);
			Arbre<String> a = arbreMap.get(s);
			
			if (time2index.containsKey(s)) {
				int i = time2index.get(s);
				m.set(i, i, 1);
				
				int j = edge2index.get(s);
				m.set (i,j,-1/r);
				
				for (Arbre<String> l:a.getChildren() ) {
					j = edge2index.get(l.getContents());
					m.set(i,j,1/r);
				}
			}
			
			int i = edge2index.get(s);
			m.set(i,i,-v/r);
		}
		
		return Math.abs(m.det());
	}
	
	// Use this to set the times and the rates
	public void setTreeState (HashMap <String,Double> times, ArrayList<String> order) {
		timeState.setTimes(times);
		timeState.setOrder(order);
		
		HashMap <String, Double> ratesMap = new HashMap <String, Double> ();
		for (String s: branchLengths.keySet()) {
			
			Arbre<String> as = arbreMap.get(s);
			if (as.isRoot())
				continue;
			String parent = as.getParent().getContents();
			double t1  = timeState.getTime(s);
			double t2  = timeState.getTime(parent);
			double v = branchLengths.get(s);
			double rate = (t2-t1)/v;
			ratesMap.put(s, rate);
		}
		rateState.setRates(ratesMap);
		
	}
	
	
	
	public void readNewick (String treeFile, Set<Taxon> inGroup)  {
		cladesMap = new HashMap <Set<String>, String> ();
		arbreMap = new HashMap <String, Arbre<String>> ();
		
		try {
//			NewickParser np = new NewickParser(IOUtils.openIn (treeFile) );
			UnrootedTree nonClockTree = UnrootedTree.fromNewick(new File(treeFile));	
			Set<Taxon> outGroup = nonClockTree.leavesSet();
			outGroup.removeAll(inGroup);
			
			// Root at an arbitrary outgroup.
			// Ensures that the ingroup forms a rooted clade.
			Arbre<Taxon> ar =  Arbre.tree2Arbre(nonClockTree.toTree (outGroup.iterator().next()));
			// Root at the parent of the LCA of the ingroup. 
			Arbre<Taxon> truncatedAr = Arbre.lowestCommonAncestor(ar, inGroup).copy();
			
			Taxon outgroup = new Taxon ("outgroup");
			Taxon root = new Taxon ("root");
			Arbre<Taxon> newAr =  Arbre.arbre(new Taxon("root"));
			newAr.addLeaves (truncatedAr.copy());
			newAr.addLeaves (Arbre.arbre(outgroup));

			Arbre<String> arOfString = Arbre.arbreToArbreOfStrings(newAr);
			this.ar = arOfString;
			
			//Tree<String> tree= np.parse();
			//ar = Arbre.tree2Arbre(tree);
			//Arbre<String> r =  ar.root ();
			
			
			Map <Arbre<String>, Set<String>> leavesMap = Arbre.leavesMap(arOfString);
			
			for (Arbre<String> a:leavesMap.keySet()){
				String label = a.getContents();
				cladesMap.put(leavesMap.get(a),label);
				arbreMap.put(label, a);
			}
			
			Map<Taxon, Double> tmpbranchLengths = new HashMap <Taxon, Double> ();
			for (Arbre<Taxon> l:truncatedAr.nodes()) {
				if (!l.isRoot() ){
					Arbre<Taxon> parent = l.getParent();
					try {
					double branchLength = nonClockTree.branchLength(l.getContents(), parent.getContents());		

					tmpbranchLengths.put (l.getContents(), branchLength);
	         } catch (Exception e)
	          {
	            System.out.println(l.getContents() + " " + parent.getContents());
	            throw new RuntimeException(e);
	          }
				} 
			}
			tmpbranchLengths.put(truncatedAr.root().getContents(),1.0);
			tmpbranchLengths.put(outgroup, 1.0);
			
			branchLengths = new HashMap<String,Double> ();
			for (Taxon l:tmpbranchLengths.keySet()) {
				branchLengths.put(l.toString(), tmpbranchLengths.get(l));
			}
	
			System.out.println(newAr.deepToString());
			
		} catch (Exception e) {
			e.printStackTrace();
		}	
	}
	
	private Set<Taxon> readInGroup (String inGroupFile) {
		HashSet <Taxon> keys = new HashSet<Taxon> ();
		for (String line:IO.i(inGroupFile)){
			StringTokenizer tok = new StringTokenizer (line,",");			
			for (;tok.hasMoreTokens(); ) {
				keys.add(new Taxon (tok.nextToken()));
			}
			break;			
		}
		return keys;
	}
 	
	public TreeState (String treeFile, String timesFile, String ratesFile, String ingroupFile) {
		Set<Taxon> inGroup  = readInGroup (ingroupFile);
		readNewick (treeFile, inGroup);
		timeState = new TimeState (ar, cladesMap, arbreMap, timesFile);
		rateState = new RateState (ar, cladesMap, arbreMap, ratesFile);
	}
	
	public TreeState (Arbre<String> ar, HashMap <Set<String>, String> cladesMap, HashMap <String, Arbre<String>> arbreMap, TimeState timeState, RateState rateState) {
		this.ar = ar;
		this.cladesMap = cladesMap;
		this.arbreMap = arbreMap;
		this.timeState = timeState;
		this.rateState = rateState;
	}
	
	public TreeState copy () {
		return new TreeState (ar, cladesMap, arbreMap, timeState, rateState);
	}
	
	public void init (Random rand)  {
		timeState.init(rand);
	}
	
	public void setTimeState (TimeState ts) {
		timeState = ts;
	}
	
	
	
	 
	public Pair<State,Double> propose (Random rand )  {
		double t = rand.nextDouble();
		TimeState currentTimeState = getTimeState ();
		TimeState proposedTimeState;
		if  (  t < TreeState.alpha ) {
			// Sample ordering and times
			proposedTimeState = currentTimeState.copy();
			proposedTimeState.setOrder(proposedTimeState.resampleOrder(rand));
			proposedTimeState.setTimes(proposedTimeState.resampleAllTimes(rand));
			
		} else {
			proposedTimeState = currentTimeState.copy();
			List<String> internalNodes = proposedTimeState.getOrder();
			int choice = rand.nextInt(internalNodes.size());
			String s = internalNodes.get(choice);
			proposedTimeState.setTimes (proposedTimeState.resampleTime(s,rand));
			// Sample only a single random time
		}
		
		double logqratio = TimeState.logTransitionProbability(proposedTimeState, currentTimeState) 
							- TimeState.logTransitionProbability(currentTimeState, proposedTimeState); 
		TreeState proposedState = this.copy();
		proposedState.setTimeState (proposedTimeState);
		
		return new Pair<State,Double>(proposedState, logqratio);
				
	}
	
	public String toString () {
		return timeState.toString();
	}
}
