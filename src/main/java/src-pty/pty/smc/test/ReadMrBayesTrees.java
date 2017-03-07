package pty.smc.test;
import java.io.*;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.regex.Pattern;

import pty.eval.SymmetricDiff;

import goblin.BayesRiskMinimizer;
import goblin.DataPrepUtils;
import goblin.Taxon;

import ma.newick.NewickParser;
import ma.newick.ParseException;
import nuts.io.IO;
import nuts.lang.StringUtils;
import nuts.util.Arbre;
import nuts.util.Counter;
import nuts.util.Tree;

public class ReadMrBayesTrees
{
  public static final Pattern selectTranslate = Pattern.compile("\\s*([0-9]+) ([^,]+)[,;]"),
                              selectTree      = Pattern.compile("\\s*tree rep[.][0-9]+ [=] (.*)$");
  public static void main(String [] args)
  {
    List<Tree<String>> trees = readTrees(args[0],10);
    Counter<Set<Set<String>>> clades = new Counter<Set<Set<String>>>();
    for (Tree<String> t : trees)
      clades.incrementCount(SymmetricDiff.clades(Arbre.tree2Arbre(t)),1.0);
    Set<Set<String>> mbr = new BayesRiskMinimizer<Set<Set<String>>>(new SymmetricDiff<Set<String>>()).findMin(clades);
    System.out.println(DataPrepUtils.newick(Arbre.arbre2Tree(SymmetricDiff.clades2arbre(mbr))));
  }
  public static List<Tree<String>> readTrees(String mrBayesOutputPath)
  {
    return readTrees(mrBayesOutputPath, 1);
  }
  public static List<Tree<String>> readTrees(String mrBayesOutputPath, int interval)
  {
    Map<String,String> translation = new HashMap<String,String>(); // int repn -> name
    List<Tree<String>> trees = new ArrayList<Tree<String>>();
    int n = 0; 
    for (String line : IO.i(mrBayesOutputPath))
      if (selectTranslate.matcher(line).matches())
      {
        List<String> matches = StringUtils.multiSelectFirstRegex(selectTranslate, line);
        translation.put(matches.get(0), matches.get(1));
      }
      else if (selectTree.matcher(line).matches()) if ((n++) % interval == 0)
      {
        String treeStr = StringUtils.selectFirstRegex(selectTree, line);
        for (String code : translation.keySet()) // hack
          for (String s : Arrays.asList("[(]", "[,]"))
            treeStr = treeStr.replaceFirst(s + code + ":" , s.charAt(1) + translation.get(code) + ":" );
//        System.out.println(treeStr);
        NewickParser np = new NewickParser(treeStr);
        try {trees.add(np.parse()); }
        catch (ParseException e) { throw new RuntimeException(e); }
      }
    return trees;
  }
}
