package pty.smc;

import fig.basic.Pair;
import goblin.Taxon;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import nuts.util.CollUtils;

public class MapLeaves {
	private final  Map<Taxon, Taxon> languageGeneMap;
	private final  Map<Taxon, Taxon> geneLanguageMap;
	
	@Override
	public String toString()
	{
	  StringBuilder result = new StringBuilder();
	  for (Taxon key : languageGeneMap.keySet())
	    result.append("" + key + "\t" + languageGeneMap.get(key) + "\n");
	  return result.toString();
	}
	
	public Set<Taxon> mapItemsNotInData(Set<Taxon> ... data)
	{
	  Set<Taxon> _data = new HashSet<Taxon>();
	  for (Set<Taxon> cur : data)
	    _data.addAll(cur);
	  _data.removeAll(languageGeneMap.keySet());
	  _data.removeAll(geneLanguageMap.keySet());
	  return _data;
	}
	public Set<Taxon> dataNotInMapItems(Set<Taxon> ... data)
	{
	  Set<Taxon> keys = CollUtils.union(languageGeneMap.keySet(), geneLanguageMap.keySet());
	  for (Set<Taxon> cur : data)
	    keys.removeAll(cur);
	  return keys;
	}
	
	public static MapLeaves parse(String filename) {
		BufferedReader br;
	    String s;
	    Map<Taxon, Taxon> languageGeneMap = new HashMap <Taxon, Taxon>();
	    Map<Taxon, Taxon> geneLanguageMap = new HashMap <Taxon, Taxon>();
	    try{
	    	br=new BufferedReader(new FileReader(filename));
	    	int linenum = 0;
	    	while((s=br.readLine())!=null){
	    		StringTokenizer tok=new StringTokenizer(s," ");
	    		String[] tmp=new String[2];
	            for(int i=0;tok.hasMoreTokens() && i<2;)
	                tmp[i++]=tok.nextToken();
	            languageGeneMap.put(new Taxon(tmp[0]), new Taxon(tmp[1]));
	            geneLanguageMap.put(new Taxon(tmp[1]), new Taxon(tmp[0]));
	                
	            linenum ++;
	    	}
	    } catch (IOException e) {
	    	e.printStackTrace();
	    }
	  return new MapLeaves(languageGeneMap, geneLanguageMap);
	}
	
	public MapLeaves(Map<Taxon, Taxon> languageGeneMap, Map<Taxon, Taxon> geneLanguageMap)
	{
	  this.languageGeneMap = CollUtils.archive(languageGeneMap);
	  this.geneLanguageMap = CollUtils.archive(geneLanguageMap);
	}
	
	public static <S,T> Set<T> mapSet(Set<S> s, Map<S,T> map) {
		HashSet<T> translatedLanguages = new HashSet<T>();
		for (S l:s) {
			T value = map.get(l);
			if (value != null)
				translatedLanguages.add (value);
		}
		return translatedLanguages;
	}
	
	public Taxon translate(Taxon l)
	{
	  Taxon 
	    t1 = languageGeneMap.get(l),
	    t2 = geneLanguageMap.get(l);
	  if (t1 != null && t2 != null)
	    throw new RuntimeException();
	  if (t1 != null) return t1;
	  if (t2 != null) return t2;
	  return null;
	}
	
	public Set<Taxon> translate (Set<Taxon> s) {
		Set<Taxon> s1 = mapSet (s, geneLanguageMap),
		              s2 = mapSet (s, languageGeneMap);
		return atMostOneEmptySet(s1,s2);
	}
	
	private static <S> Set<S> atMostOneEmptySet(Set<S> s1, Set<S> s2)
	{
	  if (!s1.isEmpty() && !s2.isEmpty()) 
	    throw new RuntimeException();
    else if (!s1.isEmpty())             return s1;
    else if (!s2.isEmpty())             return s2;
    else                                return s1; // both empty
	}
	
	public Set<Taxon> restrict(Set<Taxon> s)
	{
	  HashSet<Taxon> s1 = new HashSet<Taxon>(s),
	                    s2 = new HashSet<Taxon>(s);
	  s1.retainAll(geneLanguageMap.keySet());
	  s2.retainAll(languageGeneMap.keySet());
	  return atMostOneEmptySet(s1,s2);
	}
	/**
	 * Should use mapPairOfClades
	 * @param clades
	 * @return
	 */
  public Set<Set<Taxon>> mapClades(Set<Set<Taxon>> clades)
  {
    Set<Set<Taxon>> result = new HashSet<Set<Taxon>>();
    Set<Taxon> currentMapped;
    for (Set<Taxon> clade : clades)
      if (!(currentMapped = translate(clade)).isEmpty())
        result.add(currentMapped);
    return result;
  }
//  /**
//   * Reduces both to the domain,image of the map
//   * cleaning up empty clades, etc
//   * @param clades1
//   * @param clades2
//   * @return
//   */
//  public Pair<Set<Set<Language>>,Set<Set<Language>>> mapPairOfClades(
//      Set<Set<Language>> clades1, 
//      Set<Set<Language>> clades2)
//  {
//    Set<Set<Language>> result 
//    return null;
//  }
//  public static <S> Set<S> filterClade(Set<S> clade, Set<S> img)
//  {
//    Set<S> result = new HashSet<S>(clade);
//    result.retainAll(img);
//    return result;
//  }
  public Set<Set<Taxon>> filterClades(Set<Set<Taxon>> clades)
  {
    Set<Set<Taxon>> result = new HashSet<Set<Taxon>>();
    Set<Taxon> currentFiltered;
    for (Set<Taxon> clade : clades)
      if (!(currentFiltered = restrict(clade)).isEmpty())
        result.add(currentFiltered);
    return result;
  }

  public Map<Taxon, Taxon> getLanguageGeneMap()
  {
    return languageGeneMap;
  }

  public Map<Taxon, Taxon> getGeneLanguageMap()
  {
    return geneLanguageMap;
  }
	
//	public boolean containsAll(Set<Language> s1, Set<Language> s2) {
//		return s1.containsAll(translate(s2));				
//	}
//	
//	public boolean intersects (Set<Language> s1, Set<Language> s2) {
//		return CollUtils.intersects(s1, translate(s2));		
//	}
//	
//	public boolean contains (Set<Set<Language>> s, Set<Language> s2 ) {
//		Set<Language> newSet = translate(s2);
//		for (Set<Language> s1:s) {
//			if (s1.containsAll(newSet))
//				return true;
//		}
//		return false;
//	}
	
	
}
