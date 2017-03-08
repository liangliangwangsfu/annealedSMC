package times;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import nuts.util.Arbre;

import file.FileUtils;

public class RateState {
	private Arbre<String> ar;
	private HashMap<Set<String>, String > cladesMap;
	private HashMap <String, Arbre<String> > arbreMap;
	
	private HashMap<String, PriorType> ratesFlagMap;
	private HashMap<String, Double> ratesSDMap;
	private HashMap<String, Double> ratesMeanMap;
	private HashMap<String, Double> ratesMap;
	
	private double mean, variance;
	
	public static enum PriorType  {
		FLAT, DELTA, MOMENT
	}
			
	
	public RateState (Arbre<String> ar, HashMap <Set<String>, String> cladesMap, HashMap <String, Arbre<String>> arbreMap, String ratesFile) {
		this.ar = ar;
		this.cladesMap = cladesMap;
		this.arbreMap = arbreMap;
		this.mean = 0;
		this.variance = 1;
		readRates (ratesFile);
	}
	
	public void setMean (double mean) {
		this.mean = mean;
	}
	
	public void setVariance (double variance) {
		this.variance = variance;
	}
	
	public void setRates (HashMap <String, Double> ratesMap) {
		this.ratesMap = ratesMap;
	}
	
	public double getRate(String s) {
		return ratesMap.get(s);
	}
	
	public double logInformativePrior (double t, String s) {
		if  (ratesFlagMap.get(s) == PriorType.DELTA) {
			double u = ratesMeanMap.get (s);
			return Math.abs (t-u) <= Double.MIN_VALUE ? 0:Double.NEGATIVE_INFINITY;
		} else if (ratesFlagMap.get(s) == PriorType.FLAT) {
			return -0.5 * Math.pow(( Math.log(t) - mean) ,2)/variance - Math.log(2*Math.PI*variance);
		} else {
			double u = ratesMeanMap.get (s);
			double v=  ratesSDMap.get(s);
			v = v*v;
			return -0.5 * Math.pow(( Math.log(t)  - u) ,2)/v - Math.log(2*Math.PI*v);
		}
		
	}
	
	public double getLogDensity ()  {
		List<Arbre<String>> nodes = ar.nodes();
		double logd = 0;
		
		for (Arbre<String> s:nodes ) {
			String t = s.getContents();
			double r = ratesMap.get(t);
			
			if (!s.isRoot()) {
				double x = logInformativePrior (r,t);
				if (Double.isInfinite(x))
					return Double.NEGATIVE_INFINITY;
				else 
					logd += x;
			}
		}
		return logd;
	}

	public void readRates (String ratesFile) { 
		ratesFlagMap = new HashMap <String, PriorType> ();
		ratesMeanMap = new HashMap <String, Double> ();
		ratesSDMap = new HashMap <String, Double> ();
		ratesMap = new HashMap < String, Double> ();
		
		try { 
			HashMap <String, String> tmpMap = FileUtils.readMap(ratesFile);
			
			for (String s : tmpMap.keySet()) { 
				StringTokenizer tok = new StringTokenizer (s,",");
				HashSet <String> keys = new HashSet<String> ();
				for (;tok.hasMoreTokens(); ) {
					keys.add(tok.nextToken());
				}
				String internalLabel = cladesMap.get(keys);
				
				tok = new StringTokenizer (tmpMap.get(s),",");
				if (tok.countTokens()==1) {
					ratesFlagMap.put (internalLabel,PriorType.DELTA);
					double val = Double.parseDouble (tok.nextToken());
					ratesMeanMap.put (internalLabel, val);
					ratesSDMap.put (internalLabel, 0.0);
				} else if (tok.countTokens()==2) {
					ratesFlagMap.put (internalLabel,PriorType.MOMENT);
					double val1 = Double.parseDouble (tok.nextToken());
					double val2 = Double.parseDouble (tok.nextToken());
					ratesMeanMap.put (internalLabel, val2);
					ratesSDMap.put (internalLabel, val1);

				}				

			}
			
		} catch (Exception e)  {
			System.err.println ("Could not read from file " + ratesFile);
		}
	}

}


