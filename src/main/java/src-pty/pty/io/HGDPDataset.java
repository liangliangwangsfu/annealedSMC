/**
 * 
 */
package pty.io;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nuts.io.IO;
import nuts.lang.StringUtils;
import nuts.util.Counter;
import nuts.util.MathUtils;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import goblin.DataPrepUtils;
import goblin.Taxon;

public class HGDPDataset implements Dataset
{
  @Option public static String path = "data/hgdp/hgdp.ie.phylip";
  @Option public static int maxIndividualPerPopulation = 1;
  @Option public static ArrayList<String> restrictToPopulations = new ArrayList<String>();
  @Option public static boolean usePhylipNameTruncations = false;
  @Option public static int maxNSites = Integer.MAX_VALUE;
  
  @Option public static boolean filterSNPs = false;
  @Option public static double fractionToFilterOut = 0.9;
  
  Map<Taxon, double[][]> obs = new HashMap<Taxon,double[][]>();
  Map<Taxon, String> clu = new HashMap<Taxon, String>();
  
  Map<Taxon, Integer> popsizes = new HashMap<Taxon, Integer>();
  
  private int nChars = 2;
  private int nSites = -1;
  
  public static final Pattern POPULATION_PATTERN = Pattern.compile("^([^\\s]*)\\s+.*");
  public static final Pattern SITE_PATTERN = Pattern.compile("^[^\\s]*\\s+([0-9].*)$");
  public static final Pattern PHYLIP_CONTENT_LINE_PATTERN = Pattern.compile("^[^\\s]*\\s+[0-9].*$");
  public static final Pattern HEADERLINE = Pattern.compile("^([0-9]+\\s*)*$");
  
  public HGDPDataset()
  {
    int cIndInPop = 0;
    String prevPop = null;
    Set<String> allPops = new HashSet<String>();
    for (String line : IO.i(path))
      if (HEADERLINE.matcher(line).matches()) 
        ;
      else if (PHYLIP_CONTENT_LINE_PATTERN.matcher(line).matches())
      {
        String pop = StringUtils.selectFirstRegex(POPULATION_PATTERN, line);
        allPops.add(pop);
        if (!pop.equals(prevPop))
          cIndInPop = 0;
        String codeline = StringUtils.selectFirstRegex(SITE_PATTERN, line);
        String [] codes = codeline.split("\\s+");
        final int curNSites = Math.min(codes.length,maxNSites);
        double [][] codeArray = new double[curNSites][2];
        if (nSites == -1)
          nSites = curNSites;
        else if (nSites != curNSites)
          throw new RuntimeException();
        for (int i =0; i < nSites; i++)
        {
        	try {
          codeArray[i][0] = Double.parseDouble(codes[i]);
          codeArray[i][1] = 1.0-Double.parseDouble(codes[i]);
        	} catch (Exception e) {
        		System.out.println ("Code = " + pop + "," +  i + "," + codes[i] 
        		                     + "Prefix = " + codes[i-2] + "\t" + codes[i-1]);
        	}
        }
        Taxon lang = language(pop,cIndInPop); 
        if (cIndInPop + 1 <= maxIndividualPerPopulation && popIncluded(pop))
        {
          obs.put(lang,codeArray);
          clu.put(lang,pop); 
        }
        cIndInPop++;
        prevPop = pop;
      }
    checkPops(allPops);
    if (filterSNPs)
      obs = filterSNPs(obs);
  }
  
  private Map<Taxon, double[][]> filterSNPs(Map<Taxon, double[][]> obs)
  {
    // decide which sites to keep
    List<Integer> indicesToKeep = new ArrayList<Integer>();
    readPopulationSizes("data/hgdp/popSizes.txt");
    Counter<Integer> c = new Counter<Integer>();
    for (int s = 0; s < nSites(); s++)
      c.setCount(s, anovaTestStatistics(s));
    int requiredSites = (int) ((1.0 - fractionToFilterOut) * nSites()); // nSites() == n site before prune
    if (requiredSites < 1 || requiredSites >= nSites)
      throw new RuntimeException();
    int cSite = 0;
    for (int s : c)
    {
      indicesToKeep.add(s);
      cSite++;
      if (cSite >= requiredSites)
        break;
    }
    // construct new obs
    Map<Taxon, double[][]> newObs = new HashMap<Taxon,double[][]>();
    for (Taxon lang : obs.keySet())
    {
      final double[][] 
         cur = new double[requiredSites][],
         model=obs.get(lang);
      for (int i = 0; i < requiredSites; i++)
        cur[i] = model[indicesToKeep.get(i)];
      newObs.put(lang, cur);
    }
    // set nSites
    nSites = requiredSites;
    return newObs;
  }
  
  /*
  public static void main(String [] args) throws IOException
  {
    int n = 50000;
    HGDPDataset.path = "data/generatedExperiment/simulated.txt";
    HGDPDataset.usePhylipNameTruncations= true;
    HGDPDataset data = new HGDPDataset();
    Map<Language,double[][]> obs = data.observations();
    PrintWriter out = IOUtils.openOut("infile");
    out.println("  " + obs.size() + " " + n);
    for (int i =0 ; i < n; i++)
      out.print("2 ");
    out.println();
    for (Language lang : obs.keySet())
    {
      out.print(lang+"       ");
      double[][] cur = obs.get(lang);
      for (int i = 0; i < n; i++)
        out.print(cur[i][0]+" ");
      out.println();
    }
    out.close();
  }*/
  
  private Taxon language(String pop, int index)
  {
    if (usePhylipNameTruncations)
      pop = WalsAnn.cleanForPhylip(pop).replaceAll(" ", "");
    return new Taxon((maxIndividualPerPopulation == 1 ? pop : pop+"-"+index));
  }
  
  private void checkPops(Set<String> allPops)
  {
    for (String pop : restrictToPopulations)
      if (!allPops.contains(pop))
        LogInfo.warning("Population not recognized:" + pop);
  }

  private boolean popIncluded(String lang)
  {
    if (restrictToPopulations.size() == 0) return true;
    return restrictToPopulations.contains(lang);
  }

  public Map<Taxon, String> getReferenceClusters()
  {
    return clu;
  }

  public boolean hasReferenceClusters()
  {
    return true;
  }

  public Map<Taxon, double[][]> observations()
  {
    return obs;
  }

  public int nCharacter(int site) { return nChars; }

  public int nSites() { return nSites; }
  
  public void readPopulationSizes (String file) {
	  try {
	  BufferedReader br = new BufferedReader (new FileReader (file));
	  String s;
	  while ( (s = br.readLine() )!=null) {
		  StringTokenizer tok=new StringTokenizer(s," ");
		  popsizes.put(new Taxon(tok.nextToken()), Integer.parseInt(tok.nextToken()));			
	  
	  }
	  br.close ();
	  } catch (Exception e) {
		  e.printStackTrace ();
	  }
  }
  
  public void computeAllPairDistances () {
	Taxon[] keys = new Taxon[obs.size()];
	obs.keySet().toArray(keys);
	Counter<String> counter1 = new Counter<String>();
	Counter<String> counter2 = new Counter<String>();
	
	for (int i = 0; i < keys.length; i++)
		for (int j = 0 ; j < keys.length; j++) {
			if (i==j)
				continue;
			int n1 = popsizes.get(keys[i]);
			int n2 = popsizes.get(keys[j]);
			counter1.setCount (keys[i]+","+keys[j]+
					"=",computePairwiseDistance(obs.get(keys[i]), obs.get(keys[j])));
			
			counter2.setCount (keys[i]+","+keys[j]+
					"=",computeFST(obs.get(keys[i]), obs.get(keys[j]), n1, n2));			
	
		}

	System.out.println ("Squared distance");
	for (String s:counter1)
		System.out.println (s + counter1.getCount(s));
	
	System.out.println ("FST");
	for (String s:counter2)
		System.out.println (s + counter2.getCount(s));
	
  }
  
  public double anovaTestStatistics(int site)
  {
    double acrossPopulationsFreq = acrossPopulationsFreq(site);
    if (acrossPopulationsFreq == 0.0) return 0.0;
    double numerator = 0.0, totalPopSize = 0.0;
    for (Taxon lang : obs.keySet())
    {
      final double 
        currentPopSize = popsizes.get(lang),
        currentFreq = obs.get(lang)[site][0];
      totalPopSize += currentPopSize;
      numerator += currentPopSize * (currentFreq - acrossPopulationsFreq)
                                  * (currentFreq - acrossPopulationsFreq);
    }
    return numerator / totalPopSize / acrossPopulationsFreq;
  }
  
  private double acrossPopulationsFreq(int site)
  {
    double num = 0.0, denom = 0.0;
    for (Taxon lang : obs.keySet())
    {
      final double 
        f =  obs.get(lang)[site][0],
        n =  popsizes.get(lang);
      num += f * n;
      denom += n;
    }
    return num / denom;
  }
//  private Map<Language, Double> populationMeans(int site)
//  {
//    Map<Language, Double> result = new HashMap<Language,Double>();
//    for (Language lang : obs.keySet())
//    {
//      double num = 0.0, denom = 0.0;
//      for
//    }
//    return result;
//  }
//  
//  // main difference with computeAllPairDistance() is that it is site specific
//  // and consider all pops simultaneously
//  // used for feature (site) selection
//  private double siteSpecificFST(int site)
//  {
//    // check we know all lang's pop size
//    List<Language> miss = new ArrayList<Language>();
//    for (Language lang : obs.keySet())
//      if (!popsizes.containsKey(lang))
//        miss.add(lang);
//    if (miss.size() > 0)
//      throw new RuntimeException("Missing: " + miss);
//    final double 
//      piWithin  = siteSpecificAvgPairwiseDiffWithin (site),
//      piBetween = siteSpecificAvgPairwiseDiffBetween(site);
//    return 1.0-(piWithin+1e-10)/(piBetween+1e-10);
//  }
//  
//  private double siteSpecificAvgPairwiseDiffBetween(int site)
//  {
//    double num = 0.0, denom = 0.0;
//    final List<Language> allLangs = new ArrayList<Language>(obs.keySet());
//    for (int p = 0; p < allLangs.size(); p++)
//    {
//      final double fp = obs.get(allLangs.get(p))[site][0];
//      final int np = popsizes.get(allLangs.get(p));
//      for (int q = p+1; q < allLangs.size();q++)
//      {
//        double fq = obs.get(allLangs.get(q))[site][0];
//        final int nq = popsizes.get(allLangs.get(q));
//        denom += np*nq;
//        num += fq*nq*(1.0-fp)*np;
//        num += fp*np*(1.0-fq)*nq;
//      }
//    }
//    return num/denom;
//  }
//  private double siteSpecificAvgPairwiseDiffWithin(int site)
//  {
//    double num = 0.0, denom = 0.0;
//    for (Language lang : obs.keySet())
//    {
//      double f = obs.get(lang)[site][0];
//      int np = popsizes.get(lang);
//      denom += MathUtils.nChoose2(np);
//      num += np*np*f*(1.0-f);
//    }
//    return num/denom;
//  }
  
  
  public double computePairwiseDistance (double[][] obs1, double[][] obs2) {
	  double distance = 0;
	  
	for (int i = 0 ; i<obs1.length; i++) {
		distance += Math.pow((obs1[i][0]-obs2[i][0]), 2);
	}
	distance /= obs1.length;
	return distance;
  }
  
  public double computeFST (double[][] obs1, double[][] obs2, int n1, int n2) {
	  
	  
	  int n = n1+n2;
	  int m1 = n1*(n1-1)/2;
	  int m2 = n2*(n2-1)/2;
	  int norm;
	  norm = m1+m2;
	  double f1 = 0, f2 = 0, f = 0;
	  for (int i = 0 ; i < obs1.length; i++){
		  f1 += obs1[i][0]*(1-obs1[i][0]);
		  f2 += obs2[i][0]*(1-obs2[i][0]);
		  double p = (n1*obs1[i][0] + n2*obs2[i][0])/n;
		  f += p*(1-p);
		
	  }
	  
	  f1 = f1 * ((2.0*n1)/(n1-1)) * (((double)m1)/norm);
	  f2 = f2 * ((2.0*n2)/(n2-1)) * (((double)m2)/norm);
	  f = f  * ((2.0*n)/(n-1));
	  double fst = 1 - (f1+f2)/f;
	  
	  return fst;
	  
  }
  
  public static class PrintSiteFsts
  {
    public static void main(String [] args) 
    {
      HGDPDataset.path = "data/hgdp/frequency.pops-all.chr-all.snps-sampled.txt";
      HGDPDataset hgdp  = new HGDPDataset();
      hgdp.readPopulationSizes("data/hgdp/popSizes.txt");
      Counter<Integer> c = new Counter<Integer>();
      for (int s = 0; s < hgdp.nSites(); s++)
        c.setCount(s, hgdp.anovaTestStatistics(s));
      for (int s : c)
        System.out.println(c.getCount(s) + "\t" + hgdp.acrossPopulationsFreq(s));
    }
  }
  
  public static class PrintPruned
  {
    public static void main(String [] args)
    {
      HGDPDataset.path = "data/hgdp/frequency.pops-all.chr-all.snps-sampled.txt";
//      HGDPDataset.filterSNPs = true;
//      HGDPDataset.fractionToFilterOut = 0.9;
      HGDPDataset hgdp  = new HGDPDataset();
      Map<Taxon,double[][]> data = hgdp.observations();
      PrintWriter out = IOUtils.openOutHard("data/hgdp/forcontml.phylip");
      out.println("   " + data.keySet().size() + " " + hgdp.nSites());
      for (int i = 0; i < hgdp.nSites(); i++)
        out.print("2 ");
      out.print("\n");
      for (Taxon lang : data.keySet())
      {
        out.print(WalsAnn.cleanForPhylip(lang.toString()) + " ");
        double [][] datum = data.get(lang);
        for (int s = 0; s < datum.length; s++)
          out.print("" + datum[s][0] + " ");
        out.print("\n");
      }
      out.close();
    }
    
  }
  
  public static void main (String args[]) {
	  HGDPDataset hgdp  = new HGDPDataset();
	  hgdp.readPopulationSizes("data/hgdp/popsizes.IE.txt");
	  hgdp.computeAllPairDistances();
  }
}