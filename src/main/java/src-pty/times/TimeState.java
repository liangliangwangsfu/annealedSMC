package times;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.StringTokenizer;

import nuts.util.Arbre;
import fig.basic.Pair;
import file.FileUtils;


public class TimeState  {
	
	private Arbre<String> ar;
	private HashMap<Set<String>, String > cladesMap;
	private HashMap <String, Arbre<String> > arbreMap;
	private Map <Arbre<String>, Set<String> > nadMap;
	private Map <Arbre<String>, Set<String> > descMap;
	private Map <Arbre<String>, Set<String> > leavesMap;
	
	private List <String> order;
	
	public static enum PriorType  {
		FLAT, DELTA, BOUNDED
	}
	
	private HashMap <String, Double> timesLBMap, timesUBMap;
	private HashMap <String, PriorType> timesFlagMap;
	private HashMap <String, Double> timesMap;
	private double lambda;
	
	public List<String> getOrder () 	{
		return order;
	}
	
	public HashMap<String, Double> getTimes () {
		return timesMap;
	}

	public HashMap<String, PriorType> getTimesFlag () {
		return timesFlagMap;
	}

	public double logCalibratedPrior (double t, String s) { 
		if (timesFlagMap.get(s)==PriorType.DELTA) { 
			return t==timesLBMap.get(s)? 0 : Double.NEGATIVE_INFINITY;
		} else  if (timesFlagMap.get(s) == PriorType.BOUNDED){
			return logBoundedPrior (t, timesLBMap.get(s), timesUBMap.get(s));
		}
		return 0;
	}
	
	public double logBoundedPrior (double x, double lb, double ub)  {
		double eps = 0.05;
		
		if  (! (Double.isInfinite(lb) || Double.isInfinite(ub))) {
			double alpha = (lb/(ub-lb))*((1-eps)/(0.5*eps));
			double lambda = (1/(ub-lb)) * ((1-eps)/(0.5*eps));			
			if ( x < lb){
				return Math.log(0.5*eps) + Math.log(alpha) - Math.log(lb) 
						+ (alpha - 1) * (Math.log(x) - Math.log(lb));
			} else if (x > ub) {
				return Math.log(0.5*eps) + Math.log(lambda) - lambda * (x-ub);
			} else { 
				return Math.log (1-eps) - Math.log (ub-lb);
			}
		} else  {
			System.err.println ("Bounds must be finite");
			System.exit (1);
		}
		return 0;
	}
	
	public double sampleBoundedPrior (double lb, double ub, Random rand) {
		double eps = 0.05;
		double alpha = (lb/(ub-lb))*((1-eps)/(0.5*eps));
		double lambda = (1/(ub-lb)) * ((1-eps)/(0.5*eps));			
		
		double u  = rand.nextDouble();
		if (  u < 0.5  * eps) {
			u =  rand.nextDouble ();
			double result = Math.log(lb) + (1/alpha)*Math.log(u/(0.5*eps));
			result = Math.exp(result);
			return result;
		} else if ( u < 1  - 0.5 * eps)  {
			u =  rand.nextDouble ();
			return lb * (1-u) + ub * u;
		} else {
			u = rand.nextDouble();
			double result = -Math.log(1-(u/(0.5*eps)))/lambda + ub;
			return result;
		}
	}
	
	
	
	public double logRootPrior (double x) {
		Arbre<String> a = ar.root();
		String s = a.getContents();		
		double ub = timesUBMap.get(s);
		double lb = timesLBMap.get(s);
		
		return logBoundedPrior (x,lb,ub);
	}
	
	public double sampleRootPrior (double lb, Random rand)	{
		Arbre<String> a = ar.root();
		String s = a.getContents();
		lb = Math.max(lb, timesLBMap.get(s));
		if (!timesUBMap.containsKey(s))
			throw new RuntimeException ("Need an upper bound on the time of the root");
		
		double ub = timesUBMap.get(s);
		if (lb > ub)
			throw new RuntimeException ("Cannot sample time at root");
		double result =  rand.nextDouble() * (ub-lb) + lb ;
		return result;
	}
	
	
	public double getLogDensity() {
		// Note!!! We compute this density upto proportionality constants 
		// The constants depend on (lambda, the number of ordered trees that map to this unordered tree
		// We dont' care about these because we only want the ratio of the densities
		
		double logmarginal = 0 ;
		double logjoint =  0;
		double logcalibrated = 0;
		
		double maxtime = 0 ;
		int count =  0;
		double lb = 0;

		
		for (String s: order){
			Arbre<String> node = arbreMap.get(s);
			double t = timesMap.get (s);
			 
			if ( t  < maxtime )
				return Double.NEGATIVE_INFINITY;
			
			maxtime = t;
			
			if  ( node.isRoot ()) { 
				// SPECIFY ROOT PRIOR !!!!
				
				if (count > 0) {						
					logmarginal += count * (-lambda * lb 
								+ Math.log(1 - Math.exp(-lambda * (t-lb))) - Math.log(lambda)) ;
				}
			} else  {
				Arbre<String> parent = node.getParent();
				if ( timesMap.get(parent.getContents()) > t )
					return Double.NEGATIVE_INFINITY;
				
				if ( !timesFlagMap.containsKey(s)) {						
					logjoint += -lambda * t;
					count ++;
					continue;
				} else  {
					logjoint += -lambda * t;
					logcalibrated += logCalibratedPrior(t, s);
					
					if (count == 0) {
						logmarginal += - lambda * t;
					} else {
						logmarginal += count * (-lambda * lb 
								+ Math.log(1 - Math.exp(-lambda * (t-lb))) - Math.log(lambda)) ;
						lb = t;
						count = 0;
					}					
				}				
			}					
		}
		
		double logd =  logjoint - logmarginal + logcalibrated ;
		return logd;
	}
	

	public TimeState copy () { 
		TimeState ts = new TimeState (ar, cladesMap, arbreMap, timesMap);
		ts.setOrder(order);
		return ts;
	}

	
	public TimeState (Arbre<String> ar, HashMap <Set<String>, String> cladesMap, HashMap <String, Arbre<String>> arbreMap, HashMap<String, Double> timesMap) {
		this.ar = ar;
		this.cladesMap = cladesMap;
		this.arbreMap = arbreMap;
		this.lambda = 1;
		this.timesMap = timesMap;
		nadMap = Arbre.nadMap(ar);
		descMap = Arbre.descMap(ar);
		leavesMap = Arbre.leavesMap(ar);
		
	}
	
	public TimeState (Arbre<String> ar, HashMap <Set<String>, String> cladesMap, HashMap <String, Arbre<String>> arbreMap, String timesFile) {
		this.ar = ar;
		this.cladesMap = cladesMap;
		this.arbreMap = arbreMap;
		this.lambda = 1;
		readTimes (timesFile);
		nadMap = Arbre.nadMap(ar);
		descMap = Arbre.descMap(ar);
		leavesMap = Arbre.leavesMap(ar);
		
	}
	
	public void setLambda (double lambda) { 
		this.lambda = lambda; 
	}
	
	public void setOrder (List<String> order) {
		this.order = order;
	}
	
	
	public void setTimes (HashMap <String, Double> timesMap) {
		this.timesMap = timesMap;
	}
	
	public double getTime (String s) {
		return timesMap.get(s);
	}
	
	
	public void readTimes (String timesFile)  {
		timesFlagMap = new HashMap <String, PriorType> ();
		timesLBMap = new HashMap <String, Double> ();
		timesUBMap = new HashMap <String, Double> ();
		timesMap = new HashMap <String, Double > ();
		

		try  { 
		HashMap <String, String> tmpMap = FileUtils.readMap(timesFile);
		
		for (String s : tmpMap.keySet()) { 
			StringTokenizer tok = new StringTokenizer (s,",");
			HashSet <String> keys = new HashSet<String> ();
			for (;tok.hasMoreTokens(); ) {
				keys.add(tok.nextToken());
			}
			String internalLabel = cladesMap.get(keys);
			
			tok = new StringTokenizer (tmpMap.get(s),",");
			if (tok.countTokens()==1) {
				timesFlagMap.put (internalLabel,PriorType.DELTA);
				double val = Double.parseDouble (tok.nextToken());
				timesUBMap.put (internalLabel, val);
				timesLBMap.put (internalLabel, val);
			} else if (tok.countTokens()==2) {
				timesFlagMap.put (internalLabel,PriorType.BOUNDED);
				double val1 = Double.parseDouble (tok.nextToken());
				double val2 = Double.parseDouble (tok.nextToken());
				timesUBMap.put (internalLabel, val2);
				timesLBMap.put (internalLabel, val1);
			}				
		}
		
		if ( timesFlagMap.get(ar.root().getContents()) != PriorType.BOUNDED) {
			System.err.println ("Provide bounds on the age of the root");
				throw new RuntimeException("Provide bounds on the age of the root");
		}
		
		} catch (Exception e) {
			System.err.println ("Could not read from file " + timesFile);
			e.printStackTrace();
		}
	}
	
	public void init (Random rand) { 
		List<String> newOrder =  randomOrder (rand);
		setOrder(newOrder);
		resampleAllTimes(rand);
	}
	
	public List<String> randomOrder (Random rand) {
		List<String> result = new ArrayList<String> ();
		
		List<Arbre<String>> tmp = new ArrayList<Arbre<String>> ();
		tmp.add(ar.root());
		while (tmp.size()>0){
			int choice = rand.nextInt (tmp.size());
			Arbre<String> chosen = tmp.get(choice);
			result.add(0,chosen.getContents());
			tmp.remove(choice);
			for (Arbre<String> a:chosen.getChildren()) {
				if (!a.isLeaf())
					tmp.add(a);
			}							
		}
		
		return result;		
	}
	
	private List<String> swapOrder ( int i, int j ) { 
		List<String> result = new ArrayList<String> ();
		List<String> result1 = new ArrayList<String> ();
		
		result.addAll (order.subList(0, i));
		Arbre<String> aj = arbreMap.get(order.get(j));
		for (int k = i ; k < j; k++) { 
			String s = order.get(k);
			Arbre<String> a = arbreMap.get(s);
			if (Arbre.isAncestorOf(aj, a) )
				result.add(s);
			else
				result1.add(s);
		}
		result.add(order.get(j));
		result.addAll(result1);
		result.add(order.get(i));
		result.addAll(order.subList (j,order.size()));
		
		return result;
	}
	
	private Pair<Integer, Integer> pickRandomNadPair (Random rand) {
		int i  = rand.nextInt(order.size());
		Arbre<String> a = arbreMap.get (order.get(i));
		
		
		Set<String> nads =  nadMap.get(a);
		String[] nadsArray = new String[nads.size()];
		nads.toArray (nadsArray);
		String select = nadsArray[rand.nextInt (nadsArray.length)];
		
		int j = 0 ;
		for (String s : order) { 
			if (select.equals(s))
				break;
			j++;
		}
		return Pair.makePair(i, j);
	}
	
	public List<String> resampleOrder (Random rand) {
		Pair<Integer, Integer>	p = pickRandomNadPair (rand);
		List<String> newOrder = swapOrder (p.getFirst(), p.getSecond());
		return newOrder;
	}
	
	
	// Computes the transition probability
	// modulo constants that do not affect the MH ratio
	public static double logTransitionProbability (TimeState oldState, TimeState newState) {
		List<String> order1 = oldState.getOrder ();
		List<String> order2 = newState.getOrder();
		HashMap<String, Double> time1 = oldState.getTimes();
		HashMap<String, Double> time2 = newState.getTimes ();
		
		boolean flag = true;
		if ( order1.size() == order2.size()) {
			for (int i  = 0 ; i < order1.size(); i++)
				if (!order1.get(i).equals (order2.get(i))) {
					flag = false;
					break;
				}
					
		} else 
			flag = false;
		
		double logdensity = 0;
		if (flag)  {
			// Check to see that the times differ in exactly 
			// one coordinate
			int ind =  -1;
			for (int i = 0 ; i < order1.size(); i++) {
				double t1 = time1.get (order1.get(i));
				double t2 = time2.get (order2.get(i));
				if (Math.abs(t1-t2) > Double.MIN_VALUE ) {
					if (ind == -1)
						ind = i;
					else
						throw new RuntimeException ("Found a change in more than one time");
				}					
			}
			if (ind == -1)
				throw new RuntimeException("Times are identical");
			String s = order1.get(ind);
			logdensity = newState.logTransitionProbability(s,time2.get(s));
			
		} else {
			logdensity = newState.logTransitionProbability(time2);		
		}
		
		return logdensity;
	}
	
	
	public double logTransitionProbability (String s, double time) {
		
		for (int i = 0; i < order.size(); i ++) {
			String t =  order.get(i);
			if (t.equals(s)) {
				if (timesFlagMap.get(t) == PriorType.DELTA){
					if ( Math.abs(timesMap.get(t)-time) <Double.MIN_VALUE)
						return 0;
					else
						return Double.NEGATIVE_INFINITY;
				} else {
					
					if (i==order.size()-1) {
						return logRootPrior (time);
					} else {
						double ub = timesMap.get(order.get(i+1));
						double lb = 0;
						if (i>0)
							lb = timesMap.get(order.get(i-1));
						return Math.log (1/ub-lb);
					}
				}
				
			}
		}
		throw new RuntimeException ();
	}
	
	
	public double logTransitionProbability (HashMap <String, Double> newTimesMap) {
		double ub = Double.MAX_VALUE;
		double logdensity = 0 ;
		for (int i = order.size()-1 ;  i >= 0 ; i--) {
			String s = order.get(i);
			if (timesFlagMap.get(s) == PriorType.DELTA) {
				if ( Math.abs (newTimesMap.get(s)  - timesLBMap.get(s)) > Double.MIN_VALUE ) {
					logdensity = Double.NEGATIVE_INFINITY;
					break;
				}
			}
			if (i==order.size()-1) {
				logdensity += logRootPrior(newTimesMap.get(s));
				ub = newTimesMap.get(s);
			} else {
				double lb  = 0;
				if (i > 0)
					lb = newTimesMap.get(order.get(i-1));
				logdensity += Math.log(1/(ub-lb));
				ub = newTimesMap.get(s);
			}
		}	
		return logdensity;
	}
	
	public HashMap<String, Double> resampleAllTimes (Random rand) { 
		HashMap <String, Double> result = new HashMap <String, Double> ();
		Set<String> leaves = leavesMap.get(ar.root());
		for (String l:leaves)
			result.put (l,0.0);
	
		double ub = Double.MAX_VALUE;
		for (int i = order.size() - 1; i>=0; i--)  {
			String s = order.get(i);
			Set<String> desc = descMap.get(arbreMap.get (s));
			double min = 0 ;
			for (String t:desc )  {
				if (timesFlagMap.get(t) == PriorType.DELTA)
					min =  Math.max(timesLBMap.get(t), min) ;
			}
			
			if (timesFlagMap.get(s) == PriorType.DELTA) {
				double tmp = timesLBMap.get(s);
				
				if (tmp < min) {
					throw new RuntimeException ("Inconsistent specification of times");
				} else {
					result.put (s, tmp);
				}
					
			} else {
				if ( arbreMap.get(s).isRoot()) {
					double tmp = sampleRootPrior(min,rand);
					result.put (s,tmp);
					ub = tmp;
				} else {					
					if (min > ub)
						throw new RuntimeException ("Cannot sample times for node "  + s);
					
					double tmp = rand.nextDouble() * (ub - min)	+ min;
					result.put(s, tmp);
					ub = tmp;
				}
			}
				
		}
		
		return result;
	}
	
	
	public HashMap<String,Double> resampleTime (String s, Random rand) {
		HashMap <String, Double> result = new HashMap <String, Double> ();
		for (String key:timesMap.keySet())
			result.put(key, timesMap.get(key));
		
		if (timesFlagMap.get(s)==PriorType.DELTA)
			return result;
		
		double lb = 0;
		String root = ar.root().getContents();
		if (!timesUBMap.containsKey(root))
			throw new RuntimeException ("Need an upper bound on the time of the root");
		
		double ub = timesUBMap.get(root);
		
		int index = - 1;
		for ( int i = 0 ; i < order.size(); i++) {
			if (order.get(i) == s) {
				index=  i;
				break;
			}
		}
		
		if ( index == -1)
			throw new RuntimeException("String " + s + "not found in tree");
		
		if (index >  0)
			lb = timesMap.get(order.get(index-1));
		if (index < order.size()-1)
			ub = timesMap.get(order.get(index+1));
		
		
		double newTime = lb + (ub - lb) * rand.nextDouble();
		result.put(s, newTime);
		return result;
		
	}
	
	public String toString ()  {
		String result =  "";
		for (String s : timesMap.keySet()) {
			result += s + "\t" + timesMap.get(s);
		}
		return result;
	}
}

