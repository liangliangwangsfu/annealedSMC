package pty.io;
import java.io.*;
import java.util.*;

import fig.basic.Option;
import fig.exec.Execution;

import nuts.io.IO;
import nuts.util.EasyFormat;

/**
 * A utility to pre-process and convert Wals data
 * 
 * Can convert into PHYLIP or distance matrix (latter can be bootstrapped over sites)
 * 
 * in the latter case, distance b/w two languages is computed as:
 * 
 *  # of sites known in both lang. and with the same character values
 * ------------------------------------------------------------------
 *                # of sites known in both lang.
 * 
 * or = 1 if the 2 langs have no common known characters
 * 
 * @author bouchard
 * @deprecated use WalsDataset instead
 *
 */
public class WalsData implements Runnable
{
  @Option public String walsPath = "data/wals_data";
  @Option(gloss="Restrict to family (empty for no restriction)") public String family = "";
  @Option(gloss="If a site has more than maxNCharacters characters, reject it") public int maxNCharacters = Integer.MAX_VALUE;
  @Option(gloss="If a language has less than minNFeatures features, reject it") public int minNFeatures = 20;
  @Option public ArrayList<Operation> operations = new ArrayList<Operation>(Arrays.asList(Operation.PHYLIP));
  @Option(gloss="Seed for bootstrap") public Random rand = new Random(1); 
  @Option public int nBootstrapSamples = 100;
  
  public static WalsData parse(File dir)
  {
    WalsData result = new WalsData();
    result._parse(dir);
    return result;
  }
  
  public String getFullLanguageName(String walsCode) { return lang2Full.get(walsCode); }
  public String getFamily(String walsCode) { return lang2Family.get(walsCode); }
  public String getGenus(String walsCode) { return lang2Genus.get(walsCode); }
  
  public static enum Operation 
  { 
    /**
     * converts Wals into PHYLIP format (eg. for using with PARS)
     */
    PHYLIP {
      @Override
      public void doIt(WalsData wd)
      {
        IO.so("\n" + wd.toPhylipString(wd.maxNCharacters, wd.minNFeatures, wd.family));
      }
    }, 
    /**
     * converts Wals into a distance matrix
     */
    MTX {
      @Override
      public void doIt(WalsData wd)
      {
        System.out.println();
        wd.matrices(wd.minNFeatures, wd.nBootstrapSamples, wd.rand, wd.family);
      }
    };
    public abstract void doIt(WalsData wd);
  }
  // lang -> feature -> value
  private Map<String,Map<Integer,Integer>> data = new HashMap<String,Map<Integer,Integer>>();
  private List<Integer> features = new ArrayList<Integer>(); // one indexed!!!!
  
//  private CounterMap<Integer,Integer> featureValueCounts
  
  private List<String> langs = new ArrayList<String>();
  private Map<String, String> lang2Family = new HashMap<String,String>();
  private Map<String, String> lang2Full = new HashMap<String,String>();
  private Map<String,String> lang2Genus = new HashMap<String,String>();
  
  public void run()
  {
    _parse(new File(walsPath));
    for (Operation op : operations)
      op.doIt(this);
  }
  
  public static void main(String[] args)
  {
    Execution.run(args, new WalsData());
  }
  
  private WalsData() {}
  
  /**
   * @param root path to the wals_data directory
   */
  private void _parse(File root)
  {
    File dataPtsFile = new File(root, "datapoints.tab");
    File langFile = new File(root, "languages.tab");
//    CounterMap<String,Integer> data = new CounterMap<String,Integer>();

    parseLanguages(langFile);
    parseFeatures(dataPtsFile);
  }
  
  private void parseLanguages(File langFile)
  {
    for (String line : IO.i(langFile))
      if (!line.startsWith("wals"))
      {
        String [] fields = line.split("\\t");
        lang2Family.put(fields[0],fields[5]);
        lang2Full.put(fields[0],fields[1]);
        lang2Genus.put(fields[0],fields[4]);
      }
//    System.out.println(lang2Family);
  }


  private void parseFeatures(File dataPtsFile)
  {
    Set<Integer> _features = new HashSet<Integer>();
    for (String line : IO.i(dataPtsFile))
      if (!line.startsWith("wals_code"))
      {
        String [] fields = line.split("\\t");
        String langCode = fields[0];
        Map<Integer,Integer> type2value = new HashMap<Integer,Integer>();
        data.put(langCode, type2value);
        for (int i = 1; i < fields.length; i++)
          if (!fields[i].equals(""))
          {
            type2value.put(i, Integer.parseInt(fields[i]));
            _features.add(i);
          }
      }
    langs.addAll(data.keySet());
    Collections.sort(langs);
    features.addAll(_features);
    Collections.sort(features);
  }


  public Set<Integer> getPossibleStates(int feature)
  {
    Set<Integer> result = new HashSet<Integer>();
    for (String lang : langs)
      if (hasFeature(lang, feature))
        result.add(feature(lang,feature));
    return result;
  }
  
  public List<Integer> allFeatures() { return Collections.unmodifiableList(features); }
  
  public Set<Integer> getKnownFeatures(String lang)
  {
    Set<Integer> result = new HashSet<Integer>();
    for (int feature : allFeatures())
      if (hasFeature(lang,feature))
        result.add(feature);
    return result;
  }
  
  public double distance(String lang1, String lang2, List<Integer> features) // features might be bootstraped
  {
    List<Integer> inter = new ArrayList<Integer>();
    for (int f : features)
      if (hasFeature(lang1, f) && hasFeature(lang2,f))
        inter.add(f);
    if (inter.size() == 0) return 1.0;
    double num = 0.0;
    for (int f : inter)
      if (feature(lang1, f) != feature(lang2,f))
        num++;
    return num/inter.size();
  }
  
  public void matrices(int minNFeatures, int nBoots, Random rand,String familyRestr)
  {
//    StringBuilder result = new StringBuilder();
    List<String> langs = langWithMinNFeatures(minNFeatures, familyRestr);
    for (int i = 0; i < nBoots; i++)
    {
      System.out.print("  " + langs.size() + "\n");
      List<Integer> bsFeatures = getBootstrap(rand);
      for (String lang : langs)
      {
        System.out.print(WalsAnn.cleanForPhylip(lang));    //lang,10));
        for (String lang2 : langs)
          System.out.print(EasyFormat.fmt(distance(lang,lang2,bsFeatures))+"  ");
        System.out.print("\n");
      }
      System.out.print("\n");
    }
//    return result;
  }
  

  
  private List<Integer> getBootstrap(Random rand)
  {
    List<Integer> result = new ArrayList<Integer>();
    for (int i = 0; i < features.size(); i++)
      result.add(features.get(rand.nextInt(features.size())));
    return result;
  }


  public List<String> langWithMinNFeatures(int minNFeatures,String familyRestr)
  {
    List<String> result = new ArrayList<String>();
    for (String lang : langs)
      if (familyRestr == null || familyRestr.equals("") || lang2Family.get(lang).equals(familyRestr))
        if (getKnownFeatures(lang).size() >= minNFeatures)
          result.add(lang);
    return result;
  }


  public String toPhylipString(int maxNCharacters, int minNFeatures, String familyRestr)
  {
    List<String> langs = langWithMinNFeatures(minNFeatures, familyRestr);
    int nLangs = 0;
    StringBuilder result = new StringBuilder();
    List<Integer> features = new ArrayList<Integer>(allFeatures());
    Iterator<Integer> i = features.iterator();
    while (i.hasNext())
      if (getPossibleStates(i.next()).size() > maxNCharacters)
        i.remove();
    for (String lang : langs)
//      if (getKnownFeatures(lang).size() >= minNFeatures)
      {
        nLangs++;
        result.append(WalsAnn.fillWithSpaces(lang,10));
        for (int f : features)
          if (hasFeature(lang,f))
            result.append(feature(lang,f));
          else
            result.append("?");
        result.append("\n");
      }
//      else
//        System.err.println(lang + " had only " + getKnownFeatures(lang).size() + " known features");
    return "" + nLangs + " " + features.size() + "\n" + result.toString();
  }
  
  @Override
  public String toString()
  {
    StringBuilder result = new StringBuilder();
    result.append("wals_code\t");
    for (int i : allFeatures())
      result.append("" + i + "\t");
    for (String lang : langs)
    {
      result.append("\n" + lang);
      for (int i : allFeatures())
        if (hasFeature(lang, i))
          result.append("\t" + feature(lang,i));
        else
          result.append("\t");
    }
    return result.toString();
  }

  public int feature(String lang, Integer i)
  {
    if (!hasFeature(lang,i))
      throw new RuntimeException();
    return data.get(lang).get(i);
  }

  public boolean hasFeature(String lang, Integer i)
  {
    return data.get(lang).containsKey(i);
  }





}
