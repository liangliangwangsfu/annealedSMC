package pty.io;
import fig.basic.Pair;
import goblin.Taxon;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

import pty.io.WalsDataset.BioCharacter;
import pty.io.WalsDataset.Site;
import nuts.lang.StringUtils;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class WalsProcessingScript
{
  public Set<String> languageFamiliesToIgnore = set();
  public Set<String> featuresToIgnore = set();
  public Map<String,List<List<Integer>>> characterTranslations = map();
  public Set<String> orderFeature = set();
  
  public static final Pattern mapCmdPattern = Pattern.compile("\\s*[>]MAP[-]CHARACTERS[(]([^)]*)[)].*");
  public static final Pattern cmdPattern = Pattern.compile("\\s*[>](.*)");
  public static final Pattern featurePattern = Pattern.compile("\\s*(\\S.*\\S)\\s*");
  
  public static enum Mode { MAPCHAR, IGFEAT, IGLAN, ORDER }
  
  public WalsProcessingScript(File file)
  {
    Mode currentMode = null;
    List<List<Integer>> currentArgs = null;
    for (String line : i(file))
    {
      line = line.replaceAll("[#].*", "");
      
      if (mapCmdPattern.matcher(line).matches())
      {
        currentMode = Mode.MAPCHAR;
        currentArgs = parseMapArgs
          (StringUtils.selectFirstRegex(mapCmdPattern, line));
      }
      else if (cmdPattern.matcher(line).matches())
      {
        String parsedCmdName = StringUtils.selectFirstRegex(cmdPattern, line);
        if (parsedCmdName.equals("IGNORE-LANGUAGE-FAMILY"))
          currentMode = Mode.IGLAN;
        else if (parsedCmdName.equals("IGNORE-FEATURE"))
          currentMode = Mode.IGFEAT;
        else if (parsedCmdName.equals("ORDER-FEATURE"))
          currentMode = Mode.ORDER;
        else
          throw new RuntimeException();
      }
      else if (featurePattern.matcher(line).matches())
      {
        String feature = StringUtils.selectFirstRegex(featurePattern, line);
        if (currentMode == Mode.MAPCHAR)
        {
          characterTranslations.put(feature, currentArgs);
        }
        else if (currentMode == Mode.IGFEAT)
        {
          featuresToIgnore.add(feature);
        }
        else if (currentMode == Mode.ORDER)
        {
          orderFeature.add(feature);
        }
        else if (currentMode== Mode.IGLAN)
        {
          languageFamiliesToIgnore.add(feature);
        }
        else
          throw new RuntimeException();
      }
      else if (line.matches("\\s*"))
        ;
      else
        throw new RuntimeException();
    }
  }

  private List<List<Integer>> parseMapArgs(String str)
  {
    String [] splitCol = str.split(";");
    List<List<Integer>> result = list();
    for (String subStr : splitCol)
    {
      List<Integer> current = list();
      for (String atom : subStr.split(","))
        current.add(Integer.parseInt(atom.replaceAll("\\s","")));
      result.add(current);
    }
    return result;
  }

  public boolean ignore(Taxon t)
  {
    // find family
    String family = WalsDataset.langDB.familyMap().get(t);
    return languageFamiliesToIgnore.contains(family);
  }
  
  public boolean ignore(Site feature)
  {
    return featuresToIgnore.contains(feature.toString());
  }
  
  public boolean shouldTranslate(Site site)
  {
    return characterTranslations.keySet().contains(site.toString());
  }
  

  public Map<Pair<Site,BioCharacter>, List<Pair<Site,BioCharacter>>> cache = map();
  
  public List<Pair<Site,BioCharacter>> translate(Pair<Site,BioCharacter> item)
  {
    if (!cache.keySet().contains(item))
    {
      Site site = item.getFirst();
      List<List<Integer>> map = characterTranslations.get(site.toString());
      
      BioCharacter bioChar = item.getSecond();
      List<Integer> translation = map.get(bioChar.index);
      
      List<Pair<Site,BioCharacter>> result = list();
      for (int newSiteIdx = 0; newSiteIdx < translation.size(); newSiteIdx++)
        result.add(Pair.makePair(new Site(site.toString()+"(" + newSiteIdx + ")"), new BioCharacter("" + translation.get(newSiteIdx))));
      
      cache.put(item, result);
    }

    return cache.get(item);
  }
  
  
}
