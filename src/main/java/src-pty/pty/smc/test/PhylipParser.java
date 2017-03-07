package pty.smc.test;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.io.*;

import pty.io.HGDPDataset;

import nuts.lang.StringUtils;
import goblin.Taxon;


public class PhylipParser {
	private double[][] observations;
	private ArrayList<Taxon> names;
	private int nleaves;
	private int nsites;
	
	public PhylipParser (String filename) {
	  List<List<Double>> obs = new ArrayList<List<Double>>();
		BufferedReader br;
	    String s;
//	    ArrayList<String> lines = new ArrayList<String>();
	    names = new ArrayList<Taxon>();
	    nsites = -1;
	    
	    try{
	    	br=new BufferedReader(new FileReader(filename));
	    	int linenum = 0;
	    	while((s=br.readLine())!=null){
	    	  List<Double> currentLine = new ArrayList<Double>();
	    	  obs.add(currentLine);
	    		if (HGDPDataset.PHYLIP_CONTENT_LINE_PATTERN.matcher(s).matches()) {
//	    			StringTokenizer tok=new StringTokenizer(s," ");
//	    			String[] tmp=new String[2];
//	    			
//            for(int i=0;tok.hasMoreTokens() && i<2;)
//            	tmp[i++]=tok.nextToken();
//            names.add(new Language(tmp[0]));
//            if (tmp[1] == null)
//              throw new RuntimeException();
	    		  String name = StringUtils.selectFirstRegex(HGDPDataset.POPULATION_PATTERN, s);
	    		  names.add(new Taxon(name));
	    		  String sites = StringUtils.selectFirstRegex(HGDPDataset.SITE_PATTERN, s);
            String [] siteValues = sites.split("\\s");
            if (nsites != -1 && nsites != siteValues.length)
            {
              System.err.println ("Inconsistent number of sites in file " + filename);
              System.exit(1);
            }
            else nsites = siteValues.length;
            for (String siteValue : siteValues)
              currentLine.add(Double.parseDouble(siteValue));
            
//            lines.add(tmp[1]);
            
//            int tmpnsites = tmp[1].length();
//            if (linenum==1)
//            	nsites = tmpnsites;
//            else if (nsites!=tmpnsites) {
//            	System.err.println ("Inconsistent number of sites in file " + filename);
//            	System.exit(1);
//            }
	    		}
	            linenum ++;
	    	}
	    } catch (IOException e) {
	    	e.printStackTrace();
	    }
	    
	    nleaves = obs.size();
	    
	    observations = new double[nleaves][nsites];
	    for (int i = 0; i < nleaves; i++)
	      for (int j = 0; j < nsites; j++)
	        observations[i][j] = obs.get(i).get(j);
	    
//	    observations = new double [nleaves][];
//	    for (int i = 0 ; i < nleaves; i++) {
//	    	observations[i] = new double[nsites];
//	    	char[] tmpchars = lines.get(i).toCharArray();
//	    	
//	    	for (int j = 0 ; j < nsites; j++) {
//	    		if ( tmpchars[j] == '?')
//	    			observations[i][j] = -1;
//	    		else
//	    			observations[i][j] = Character.getNumericValue(tmpchars[j]);
//	    	}
//	    }
	}
	
	public String toString () {
		String out = "";
		for (int i = 0 ; i < nleaves; i++) {
			for (int j = 0 ; j < nsites; j++) {
				out += observations[i][j] + "\t";
			}
			out += "\n";
		}
		
		return out;
	}
	
	public int getNumberOfLeaves () {
			return nleaves;
	}
	
	public int getNumberofSites () {
			return nsites;
	}
	
	
	public double[][] getDoubleObservations () {
		return observations;
	}
	
	public int [][] getObservations ()
	{
	  int[][] result = new int[observations.length][];
	  for (int i = 0; i < observations.length; i++)
	  {
	    result[i] = new int[observations[i].length];
	    for (int j =0 ; j < result[i].length; j++)
	      result[i][j] = (int) observations[i][j];
	  }
	  return result;
	}
	

	public ArrayList<Taxon> getLeafNames () {
		return names;
	}
	
	public static void main (String args[]) {
		PhylipParser p = new PhylipParser(args[0]);
		System.out.println (p.toString());
	}
	
}
