/**
 * 
 */
package pty.io;

import fig.basic.IOUtils;
import goblin.Taxon;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;


import ma.newick.NewickParser;
import ma.newick.ParseException;
import nuts.io.IO;
import nuts.lang.StringUtils;
import nuts.util.Arbre;
import nuts.util.Tree;
import nuts.util.Arbre.ArbreMap;

/**
 * A utility that substitute wals codes found in the standard in by
 * the full language name and corresponding family and genus and
 * print it in the std out
 * @author bouchard
 */
public class WalsAnn
{
  public static void main(String [] args) throws IOException, ParseException
  {
    
    List<String> walsPaths = Arrays.asList(
        "/home/eecs/bouchard/ptychodus/data/wals_data",
        "/Users/bouchard/w/ptychodus/data/wals_data");
    
    WalsData data = null;
    
    for (String path : walsPaths)
      try {
        data = WalsData.parse(new File(path));
      } catch (Exception e) {}
    
//    if (args.length >= 1)
//      walsPath = args[0];
//    WalsData data = WalsData.parse(new File(walsPath));
//    if (args.length < 2 || !args[1].equals("t"))
      for (String line : IO.si())
        IO.so(annotate(data,line));
//    else
//    {
      
//      StringBuilder b = new StringBuilder();
////      for (String line : IO.si())
////        b.append(line + "\n");
//      NewickParser np = new NewickParser(IOUtils.openIn(args[0]));
//      Tree<String> tree = np.parse();
//      IO.so(annotate(data, Arbre.tree2Arbre(tree).preOrderMap(new ArbreMap<String,Language>() {
//        @Override
//        public Language map(Arbre<String> currentDomainNode)
//        {
//          return new Language(currentDomainNode.getContents());
//        }
//      }), false));
//    }
  }
  private static final Pattern selector = Pattern.compile("[^a-z]([a-z]{3})[^a-z]");
  
  public static String annotate(final WalsData data, Arbre<Taxon> a, final boolean showFamily)
  {
    return a.preOrderMap(new ArbreMap<Taxon,String>() {
      @Override
      public String map(Arbre<Taxon> currentDomainNode)
      {
        if (currentDomainNode.getContents() == null) return "";
        return processWalsCode(data, currentDomainNode.getContents().toString(), showFamily);
      }
    }).deepToString();
  }
  
  public static String annotate(WalsData data, String line)
  {
    for (String match : StringUtils.selectRegex(selector, line))
    {
      String sub = data.getFullLanguageName(match); //lang2Full.get(match);
      if (sub != null)
        line = line.replaceAll(match, 
            processWalsCode(data, match, false));
      else 
        System.err.println("Warning: unknown code: " + match);
    }
    return line;
  }
  
  private static String processWalsCode(WalsData data, String walscode, boolean showFamily)
  {
    String sub = data.getFullLanguageName(walscode);
    if (sub != null)
      return
          cleanedLangName(sub, false);
//    +
//          "_" + (showFamily ? cleanedLangName(data.getFamily(walscode),false) + "/" : "")+
//          cleanedLangName(data.getGenus(walscode),false);
    else return walscode;
  }
  
  public static String cleanForPhylip(String s)
  {
    return fillWithSpaces(WalsAnn.cleanedLangName(s,true,10),10);
  }
  
  public static String fillWithSpaces(String s, int n)
  {
    if (s.length() > n) 
      return s.substring(0,n);
    int delta = n - s.length();
    for (int i = 0; i < delta; i++)
      s += " ";
    return s;
  }
  
  // clip language names to 10 character for the stupid phylip input format
  // Warning: with WALS, this will create duplicates!
  public static String cleanedLangName(String s)
  { return cleanedLangName(s, true); }
  // useful for removing parentheses and accents in names and optionally clipping to 10
  public static String cleanedLangName(String s, boolean clip)
  {
    return cleanedLangName(s, clip, 10);
  }
  public static String cleanedLangName(String s, boolean clip, int clipL)
  {
    // keep only letters
    String result = "";
    for (char c : s.toCharArray())
      if (("" + c).matches("[A-Za-z0-9]"))
        result += c;
    if (result.length() == 0)
      throw new RuntimeException("orig:" + s);
    // up to length 10 (for phylip)
    if (clip)
      return result.substring(0,Math.min(clipL,result.length())); 
    else return result;
  }
}