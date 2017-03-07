package pty.io;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import fig.exec.Execution;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import pty.smc.MapLeaves;

import nuts.io.IO;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.CounterMap;
import nuts.util.Indexer;

/**
 * Actually, more general than WALS
 * TODO: decouple better
 * @author bouchard
 *
 */
public class WalsDataset implements Dataset
{
  @Option public static String scriptPath = "/Users/bouchard/Documents/workspace/evolvere/data/wals-preprocessing-script";
  @Option public static boolean useWalsCodeForLanguages = true;
  @Option public static String walsPath = "data/wals_data";
  @Option public static String languageFamilyRestriction = "Indo-European";
  @Option public static ArrayList<String> languageListRestriction = new ArrayList<String>();
  @Option public static double sitesFractionThreshold = 0.25; 
  @Option public static double charsFractionThreshold = 0.25; 
  @Option public static int languageCountThreshold = 10;
  @Option public static boolean useFamilyAsRef = false;
  @Option public static ArrayList<WalsCorpusOperation> preprocessingSteps 
    = new ArrayList<WalsCorpusOperation>(Arrays.asList(
        WalsCorpusOperation.FAMILY_RESTRICT,
        WalsCorpusOperation.UNDERDOCUMENTED_LANGS,
        WalsCorpusOperation.UNDERUSED_SITES,
        WalsCorpusOperation.BINARIZE,
        WalsCorpusOperation.REMOVE_DEGENERATE_SITES));
  
  public static WalsDataset getPreprocessedCorpus()
  {
    WalsDataset initial = null;
//    if (useLexData)
//    {
//      LogInfo.track("Preprocessing Warnow's lexical data", true);
//      WarnowParser parser = new WarnowParser(new File(warnowPath));
//      initial = parser.getDataset();
//    }
//    else
//    {
      LogInfo.track("Preprocessing WALS", true);
      Parser parser = new Parser(new File(walsPath));
      initial = parser.getDataset();
//    }
    LogInfo.logs("Initial: " + initial.summary());
    WalsDataset result = getProcessedDataset(initial, preprocessingSteps);
    LogInfo.end_track();
    return result;
  }
  
  public Map<Taxon,double[][]> observations()
  {
    return Dataset.DatasetUtils.convert(toObservationArrays(unknownCode), this);
  }
  

  public boolean hasReferenceClusters() { return true; }
  public Map<Taxon,String> getReferenceClusters()
  {
    return (useFamilyAsRef ? langDB.familyMap() : langDB.genusMap());
  }
  
  private final Map<Pair<Taxon,Site>,BioCharacter> data;
  
  private WalsDataset(Map<Pair<Taxon,Site>,BioCharacter> data)
  {
    this.data = CollUtils.archive(data);
  }
  
  public String summary() 
  {
    return "" + allLanguages().size() + " languages and " + allSites().size() + " sites";
  }
  public String toString(Taxon lang)
  {
    StringBuilder result = new StringBuilder();
    if (!allLanguages().contains(lang)) throw new RuntimeException();
    result.append(lang + "\n");
    for (Site site : knownSites(lang))
      result.append("\t" + site + "\t" + get(lang,site) + "\n");
    return result.toString();
  }
  @Override
  public String toString()
  {
    StringBuilder result = new StringBuilder();
    result.append(summary()+"\n"); 
    for (Taxon lang : allLanguages())
      result.append(toString(lang));
    return result.toString();
  }
  
  public String toPhylip()
  {
    final int unk = -1;
    Map<Taxon,int[]> arrays = toObservationArrays(unk);
    StringBuilder result = new StringBuilder();
    result.append("" + allLanguages().size() + " " + allSites().size() + "\n");
    for (Taxon lang : arrays.keySet())
    {
      // use langcode for making sure they are unique when truncated
      String langCode = langDB.findWalsCode(lang);
      result.append(WalsAnn.cleanForPhylip(langCode));
      for (int site : arrays.get(lang))
        if (site == unk)
          result.append("?");
        else
          result.append(site);
      result.append("\n");
    }
    return result.toString();
  }
  
  private Map<Taxon,int[]> toObservationArrays(int unknownCode)
  {
    Map<Taxon,int[]> result = new HashMap<Taxon,int[]>();
    Indexer<Site> siteIndexer = siteIndexer();
    Map<Site,Indexer<BioCharacter>>  charIndexers = charIndexers();
    for (Taxon lang : allLanguages())
    {
      int [] current = new int[siteIndexer.size()];
      for (int s = 0; s < current.length; s++)
      {
        Site site = siteIndexer.i2o(s);
        if (!isKnown(lang,site))
          current[s] = unknownCode;
        else
        {
          int charIndex = charIndexers.get(site).o2i(get(lang,site));
          current[s] = charIndex;
        }
      }
      result.put(lang,current);
    }   
    return result;
  }
  
  private transient Indexer<Site> _siteIndexer = null;
  public Indexer<Site> siteIndexer()
  {
    if (_siteIndexer != null) return _siteIndexer;
    List<Site> sorted = new ArrayList<Site>(allSites());
    Collections.sort(sorted);
    _siteIndexer = new Indexer<Site>(sorted);
    return _siteIndexer;
  }
  private transient Map<Site,Indexer<BioCharacter>> _charIndexers = null;
  public Map<Site,Indexer<BioCharacter>> charIndexers()
  {
    if (_charIndexers != null) return _charIndexers;
    _charIndexers = new HashMap<Site,Indexer<BioCharacter>>();
    for (Site site : allSites())
    {
      List<BioCharacter> sorted = new ArrayList<BioCharacter>(allCharacters(site));
      try{
        Collections.sort(sorted);
      } catch (Exception e)
      {
        System.out.println(e);
      }
      _charIndexers.put(site, new Indexer<BioCharacter>(sorted));
    }
    return _charIndexers;
  }
  
  public static WalsProcessingScript getScript()
  {
    return new WalsProcessingScript(new File(scriptPath));
  }
  
  public static Set<Taxon> getFamilyRestriction() 
  { 
    if (languageFamilyRestriction == null || languageFamilyRestriction.equals(""))
      return langDB.family.keySet();
    Set<Taxon> result = new HashSet<Taxon>();
    for (Taxon lang : langDB.family.keySet())
      if (langDB.family.get(lang).equals(languageFamilyRestriction))
        result.add(lang);
    return result;
  }
  
  public static Set<Taxon> getListRestriction() 
  { 
    if (languageListRestriction.size() == 0)
      return langDB.family.keySet();
    Set<Taxon> result = new HashSet<Taxon>();
    for (String langStr : languageListRestriction)
    {
      Taxon cur = new Taxon(langStr);
      if (!langDB.family.keySet().contains(cur))
        throw new RuntimeException();
      result.add(cur);
    }
    return result;
  }
  
  public static LanguageDatabase langDB = null;
  public static void loadLanguageDatabase(File walsDir)
  {
    if (langDB != null) return;
//    if (useLexData)
//      langDB = new LanguageDatabase(new File(lexDataLangs));
//    else
      langDB = new LanguageDatabase(new File(walsDir, "languages.tab"));
  }
  public static class LanguageDatabase
  {
    private final Map<Taxon, String> genus = new HashMap<Taxon,String>();
    private final Map<Taxon, String> family = new HashMap<Taxon,String>();
    private final Map<String, Taxon> walsCode2Language = new HashMap<String,Taxon>();
    public Taxon walsCode2Language(String wals) { return walsCode2Language.get(wals); }
    public Map<Taxon, String> genusMap() { return Collections.unmodifiableMap(genus); }
    public Map<Taxon, String> familyMap() { return Collections.unmodifiableMap(family); }
    private LanguageDatabase(File languageFile)
    {
      for (String line : IO.i(languageFile))
        if (!line.matches("^wals.*"))
        {
          String [] fields = line.split("\\t");
          String walsCode = fields[0];
          Taxon lang = new Taxon(useWalsCodeForLanguages ? 
              walsCode : 
              WalsAnn.cleanedLangName(fields[1],false));
          genus.put(lang, fields[4]);
          family.put(lang, fields[5]);
          walsCode2Language.put(walsCode, lang);
        }
    }
    public String findWalsCode(Taxon lang)
    {
      for (String walsCode : walsCode2Language.keySet())
        if (walsCode2Language.get(walsCode).equals(lang))
          return walsCode;
      return null;
    }
  }
  
  public static WalsDataset languageRestrict(WalsDataset dataset, Set<Taxon> restr)
  {
    Map<Pair<Taxon,Site>,BioCharacter> newData =
      new HashMap<Pair<Taxon,Site>,BioCharacter>();
//    Iterator<Pair<Language,Site>> i = newData.keySet().iterator();
    for (Pair<Taxon,Site> key : dataset.data.keySet()) 
      if (restr.contains(key.getFirst()))
        newData.put(key, dataset.data.get(key));
    return new WalsDataset(newData);
  }
  
  public static void main(String [] args)
  {
    Execution.run(args, new WalsDatasetApp(), "wals", WalsDataset.class);
  }
  public static class WalsDatasetApp implements Runnable
  {
    public void run()
    {
      preprocessingSteps 
      = new ArrayList<WalsCorpusOperation>(Arrays.asList(WalsCorpusOperation.SCRIPT, WalsCorpusOperation.UNDERDOCUMENTED_LANGS));
      WalsDataset dataset = WalsDataset.getPreprocessedCorpus();
      IO.so(dataset.toPhylip());
      //
      MapLeaves ml = MapLeaves.parse("data/world-language-gene-map.txt");
      Set<Taxon> restr = ml.getLanguageGeneMap().keySet();
      Map<Taxon,String> labels = new HashMap<Taxon,String>(langDB.genusMap());
//      labels.keySet().retainAll(restr);
      Map<String,Set<Taxon>> inv = CollUtils.invert(labels);
      for (String str : inv.keySet())
        System.out.println("" + str + ":\t" + inv.get(str));
      System.out.println("---");
      labels = new HashMap<Taxon,String>(langDB.familyMap());
//      labels.keySet().retainAll(restr);
      inv = CollUtils.invert(labels);
      for (String str : inv.keySet())
        System.out.println("" + str + ":\t" + inv.get(str));
    }
  }
  
  public static WalsDataset getProcessedDataset(WalsDataset initial, List<WalsCorpusOperation> list)
  {
    WalsDataset current = initial;
    for (WalsCorpusOperation op : list)
    {
      current = op.apply(current);
      LogInfo.logs("After " + op.summary() + ": " + current.summary());
    }
    return current;
  }
  
  private static class Parser
  {
    private final Map<Integer,Site> sites = new HashMap<Integer,Site>();
    private final Map<Site,Map<Integer,BioCharacter>> chars = new HashMap<Site,Map<Integer,BioCharacter>>();
    private final Map<Pair<Taxon,Site>,BioCharacter> data = new HashMap<Pair<Taxon,Site>,BioCharacter>();
    public Parser(File file)
    {
      loadLanguageDatabase(file);
      parseSites(new File(file, "features.tab"));
      parseBioCharacters(new File(file, "values.tab"));
      parseDataPoints(new File(file, "datapoints.tab"));
    }
    public WalsDataset getDataset() { return new WalsDataset(data); }
    private void parseDataPoints(File file)
    {
      List<Integer> featureCodes = null;
      for (String line : IO.i(file))
        if (!line.matches("^wals_code.*"))
        {
          String [] fields = line.split("\\t");
          Taxon lang = langDB.walsCode2Language(fields[0]);
//          if (fields.length != featureCodes.size() + 1)
//            throw new RuntimeException(line + "\n" + fields.length + " vs " + (  featureCodes.size() + 1));
          for (int i = 1; i < fields.length; i++)
            if (!fields[i].equals("")) // not missing
            {
              int siteId = featureCodes.get(i);
              Site site = sites.get(siteId);
              int bcId = Integer.parseInt(fields[i]);
              Map<Integer,BioCharacter> cMap = chars.get(site); 
              BioCharacter bc = cMap.get(bcId);
              if (bc == null)
                throw new RuntimeException();
              data.put(Pair.makePair(lang,site), bc);
            }
        }
        else
        {
          featureCodes = new ArrayList<Integer>();
          String [] fields = line.split("\\t");
          featureCodes.add(null);
          for (int i = 1; i < fields.length; i++)
            featureCodes.add(Integer.parseInt(fields[i]));
        }
    }

    private void parseSites(File file)
    {
      for (String line : IO.i(file))
        if (!line.matches("^id.*"))
        {
          String [] fields = line.split("\\t");
          int site = Integer.parseInt(fields[0]);
          sites.put(site, new Site(fields[1]));
        }
    }

    private void parseBioCharacters(File file)
    {
      for (String line : IO.i(file))
        if (!line.matches("^feature_id.*"))
        {
          String [] fields = line.split("\\t");
          int _site = Integer.parseInt(fields[0]),
              bc   = Integer.parseInt(fields[1]);
          Site site = sites.get(_site);
          Map<Integer,BioCharacter> current = chars.get(site);
          if (current == null)
          {
            current = new HashMap<Integer,BioCharacter>();
            chars.put(site,current);
          }
          current.put(bc, new BioCharacter(fields[2], bc - 1));
        }
    }
  }

  public static enum WalsCorpusOperation
  {
    FAMILY_RESTRICT {
      @Override public String summary()
      {
        return super.toString() + "(" + languageFamilyRestriction + ")";
      }
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        return languageRestrict(dataset, WalsDataset.getFamilyRestriction());
      }
    },
    LIST_RESTRICT {
      @Override public String summary()
      {
        return super.toString() + "(" + languageListRestriction + ")";
      }
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        return languageRestrict(dataset, WalsDataset.getListRestriction());
      }
    },
    UNDERDOCUMENTED_LANGS {
      @Override public String summary()
      {
        return super.toString() + "(" + languageCountThreshold + ")";
      }
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        // identify good langs
        Set<Taxon> restr = new HashSet<Taxon>();
        for (Taxon lang : dataset.allLanguages())
          if (dataset.knownSites(lang).size() > languageCountThreshold)
            restr.add(lang);
        return languageRestrict(dataset, restr);
      }
    },
    UNDERUSED_CHARS {
      @Override public String summary()
      {
        return super.toString() + "(" + charsFractionThreshold + ")";
      }
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        // identify good characters
        double nLangs = dataset.allLanguages().size();
        Counter<Pair<Site,BioCharacter>> charCounter = dataset.getBioCharacterCounts();
        // keep only those
        Map<Pair<Taxon,Site>,BioCharacter> newData =
          new HashMap<Pair<Taxon,Site>,BioCharacter>();
        for (Pair<Taxon,Site> key : dataset.data.keySet()) 
        {
          BioCharacter value = dataset.data.get(key);
          Pair<Site,BioCharacter> key2 = Pair.makePair(key.getSecond(), value);
          if (charCounter.getCount(key2) / nLangs > charsFractionThreshold) 
            newData.put(key, dataset.data.get(key));
        }
        return new WalsDataset(newData);
      }
    },
    SCRIPT {

      @Override
      public WalsDataset apply(WalsDataset dataset)
      {
        Map<Pair<Taxon,Site>,BioCharacter> newData =
          new HashMap<Pair<Taxon,Site>,BioCharacter>();
        WalsProcessingScript script = WalsDataset.getScript();
        
        loop:for (Pair<Taxon,Site> key : dataset.data.keySet()) 
        {
          final Site site = key.getSecond(); 
          final Taxon taxon = key.getFirst();
          if (script.ignore(taxon)) continue loop;
          if (script.ignore(site)) continue loop;
          final BioCharacter value = dataset.data.get(key);
          if (script.shouldTranslate(site))
            for (Pair<Site,BioCharacter> translated : script.translate(Pair.makePair(site, value)))
              newData.put(Pair.makePair(taxon, translated.getFirst()), translated.getSecond());
          else
            newData.put(key, value);
        }
        return new WalsDataset(newData);
      }

      @Override
      public String summary() { return toString(); }
    },
    ENCODE_BIN {
      public final BioCharacter  
        ZERO = new BioCharacter("ZERO"),
        ONE  = new BioCharacter("ONE");
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        Map<Pair<Taxon,Site>,BioCharacter> newData =
          new HashMap<Pair<Taxon,Site>,BioCharacter>();
        for (Pair<Taxon,Site> key : dataset.data.keySet()) 
        {
          final BioCharacter value = dataset.data.get(key);
          // set all
          for (BioCharacter otherValue : dataset.allCharacters(key.getSecond()))
          {
            // new site:
            Site newSite = new Site("" + key.getSecond() + "=" + otherValue);
            Pair<Taxon,Site> newKey = Pair.makePair(key.getFirst(), newSite);
            if (value.equals(otherValue))
              newData.put(newKey, ONE);
            else
              newData.put(newKey, ZERO);
          }
        }
        return new WalsDataset(newData);
      }
      @Override
      public String summary() { return toString(); }
    },
    BINARIZE {
      public final BioCharacter SET_OF_COLLAPSED_NON_MODE_VALUES 
        = new BioCharacter("SET_OF_COLLAPSED_NON_MODE_VALUES");
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        // get cound info
        CounterMap<Site,BioCharacter> charCounter = toCounterMap(dataset.getBioCharacterCounts());
        Map<Pair<Taxon,Site>,BioCharacter> newData =
          new HashMap<Pair<Taxon,Site>,BioCharacter>();
        for (Pair<Taxon,Site> key : dataset.data.keySet()) 
        {
          BioCharacter value = dataset.data.get(key);
          if (charCounter.getCounter(key.getSecond()).argMax().equals(value))
            newData.put(key,value);
          else
            newData.put(key,SET_OF_COLLAPSED_NON_MODE_VALUES);
        }
        return new WalsDataset(newData);
      }
      @Override
      public String summary() { return toString(); }
    },
    REMOVE_DEGENERATE_SITES {
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        CounterMap<Site,BioCharacter> charCounter = toCounterMap(dataset.getBioCharacterCounts());
        Map<Pair<Taxon,Site>,BioCharacter> newData =
          new HashMap<Pair<Taxon,Site>,BioCharacter>();
        for (Pair<Taxon,Site> key : dataset.data.keySet()) 
          if (charCounter.getCounter(key.getSecond()).size() >= 2)
            newData.put(key,dataset.data.get(key));
        return new WalsDataset(newData);
      }
      @Override
      public String summary() { return toString(); }
    },
    REMOVE_SIGN {
      @Override public WalsDataset apply(WalsDataset dataset)
      {
        Map<Pair<Taxon,Site>,BioCharacter> newData =
          new HashMap<Pair<Taxon,Site>,BioCharacter>();
  //      Iterator<Pair<Language,Site>> i = newData.keySet().iterator();
        for (Pair<Taxon,Site> key : dataset.data.keySet()) 
          if (!langDB.familyMap().get(key.getFirst()).equals("other"))
            newData.put(key, dataset.data.get(key));
        return new WalsDataset(newData);
      }
      @Override
      public String summary() { return toString(); }
    },
    UNDERUSED_SITES {
      @Override public String summary()
      {
        return super.toString() + "(" + sitesFractionThreshold + ")";
      }
      @Override
      public WalsDataset apply(WalsDataset dataset)
      {
        // identify good sites
        double nLangs = dataset.allLanguages().size();
        Counter<Site> siteCounter = dataset.getSiteCounts();
        // keep only those
        Map<Pair<Taxon,Site>,BioCharacter> newData =
          new HashMap<Pair<Taxon,Site>,BioCharacter>();
        for (Pair<Taxon,Site> key : dataset.data.keySet()) 
          if (siteCounter.getCount(key.getSecond()) / nLangs > sitesFractionThreshold)             
            newData.put(key, dataset.data.get(key));
        return new WalsDataset(newData);
      }
    };
    public abstract WalsDataset apply(WalsDataset dataset);
    public abstract String summary();
  }
  
  public static <S,T> CounterMap<S,T> toCounterMap(Counter<Pair<S,T>> c)
  {
    CounterMap<S,T> result = new CounterMap<S,T>();
    for (Pair<S,T> key : c.keySet())
      result.setCount(key.getFirst(), key.getSecond(), c.getCount(key));
    return result;
  }
  
  public Counter<Pair<Site,BioCharacter>> getBioCharacterCounts()
  {
    Counter<Pair<Site,BioCharacter>> result = 
      new Counter<Pair<Site,BioCharacter>>();
    for (Pair<Taxon,Site> key : data.keySet())
      result.incrementCount(Pair.makePair(key.getSecond(),data.get(key)),1.0);
    return result;
  }
  /**
   * number of known 
   * @return
   */
  public Counter<Site> getSiteCounts()
  {
    Counter<Site> result = 
      new Counter<Site>();
    for (Pair<Taxon,Site> key : data.keySet())
      result.incrementCount(key.getSecond(),1.0);
    return result;
  }
  
  private transient Set<Taxon> _allLanguages = null;
  public Set<Taxon> allLanguages()
  {
    if (_allLanguages != null) return _allLanguages;
    _allLanguages = new HashSet<Taxon>();
    for (Pair<Taxon,Site> key : data.keySet())
      _allLanguages.add(key.getFirst());
    return _allLanguages;
  }
  
  private transient Set<Site> _allSites = null;
  public Set<Site> allSites()
  {
    if (_allSites != null) return _allSites;
    _allSites = new HashSet<Site>();
    for (Pair<Taxon,Site> key : data.keySet())
      _allSites.add(key.getSecond());
    return _allSites;
  }
  
  private transient Map<Site,Set<BioCharacter>> _allCharacters = new HashMap<Site,Set<BioCharacter>>();
  public Set<BioCharacter> allCharacters(Site site)
  {
    if (!allSites().contains(site)) throw new RuntimeException();
    if (_allCharacters.get(site) != null) return _allCharacters.get(site);
    for (Pair<Taxon,Site> key : data.keySet())
      CollUtils.getNoNullSet(_allCharacters,key.getSecond()).add(data.get(key));
    return _allCharacters.get(site);
  }
  
  public boolean isKnown(Taxon lang, Site site)
  {
    if (!allLanguages().contains(lang) ||
        !allSites().contains(site))
      throw new RuntimeException();
    return data.containsKey(Pair.makePair(lang,site));
  }
  
  private transient Map<Taxon,Set<Site>> _knownCharacter = new HashMap<Taxon,Set<Site>>();
  public Set<Site> knownSites(Taxon lang)
  {
    if (!allLanguages().contains(lang)) throw new RuntimeException();
    if (_knownCharacter.containsKey(lang)) return _knownCharacter.get(lang);
    Set<Site> result = new HashSet<Site>();
    for (Pair<Taxon,Site> key : data.keySet())
      if (key.getFirst().equals(lang))
        result.add(key.getSecond());
    _knownCharacter.put(lang,result);
    return result; 
  }
  
  public BioCharacter get(Taxon lang, Site site)
  {
    if (!isKnown(lang,site))
      throw new RuntimeException();
    return data.get(Pair.makePair(lang,site));
  }
  
  
  public static class Site implements Serializable, Comparable<Site>
  {
    private static final long serialVersionUID = 1L;
    private final String string;
    public Site(String string) 
    { 
      if (string == null || string.equals(""))
        throw new RuntimeException("Name of site should be nontrivial");
      this.string = string; 
    }
    @Override
    public boolean equals(Object obj)
    {
      if (!(obj instanceof Site)) return false;
      return this.string.equals(((Site) obj).string);
    }
    @Override
    public int hashCode()
    {
      return string.hashCode();
    }
    @Override
    public String toString()
    {
      return string;
    }
    public int compareTo(Site arg0)
    {
      return this.string.compareTo(arg0.string);
    }
  }
  public static  class BioCharacter implements Serializable, Comparable<BioCharacter>
  {
    private static final long serialVersionUID = 1L;
    private final String string;
    public final int index; // zero indexed
    public BioCharacter(String string) 
    { 
      this(string, -1);
    }
    public BioCharacter(String string, int index) 
    { 
//      if (index < 0) throw new RuntimeException();
      if (string == null || string.equals(""))
        throw new RuntimeException("Name of character should be nontrivial");
      this.string = string; 
      this.index = index;
    }
    @Override
    public boolean equals(Object obj)
    {
      if (!(obj instanceof BioCharacter)) return false;
      return this.string.equals(((BioCharacter) obj).string);
    }
    @Override
    public int hashCode()
    {
      return string.hashCode();
    }
    @Override
    public String toString()
    {
      return string;
    }
    public int compareTo(BioCharacter arg0)
    {
      return this.string.compareTo(arg0.string);
    }
    
  }

  public int nCharacter(int site)
  {
    final int result = allCharacters(siteIndexer().i2o(site)).size();
    return result;
  }

  public int nSites()
  {
   return allSites().size();
  }
}
