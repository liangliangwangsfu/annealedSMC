package pty.io;
import java.io.*;
import java.util.*;

import nuts.io.IO;
import nuts.util.Indexer;

public class ConvertWarnowData
{
  public static void main(String [] args)
  {
    StringBuilder values = new StringBuilder();
    String [] langNames = null;
    List<List<String>> entries = new ArrayList<List<String>>();
    int lineNumber = 0;
    System.out.println("features.tab:");
    for (String line : IO.i("data/lexical/cphl/withoutTypo.txt"))
    {
      String [] fields = line.split(" ");
      if (lineNumber == 0)
        langNames = fields;
      else
      {
        List<String> current = new ArrayList<String>();
        entries.add(current);
        for (int i = 0; i < fields.length; i++)
          if (i > 1)
            current.add("" + (Integer.parseInt(fields[i])));
        int featureId = (lineNumber);
        System.out.println("" + featureId + "\t" + fields[1]);
        Set<String> cvals = new HashSet<String>();
        for (int i = 2; i < fields.length; i++)
          cvals.add(fields[i]);
        for (String j : cvals)
          values.append("" + featureId + "\t" + j + "\tf"+featureId+"v"+j+"\tf"+featureId+"v"+j+"\n");
      }
      lineNumber++;
    }
    System.out.println("---");
    System.out.println("values.tab:");
    System.out.println(values);
    System.out.println("---");
    System.out.println("datapoints.tab:");
    System.out.print("wals_code\t");
    for (int i = 1; i <= 421; i++)
      System.out.print("" + i+"\t");
    System.out.println();
    for (int l = 0; l < langNames.length; l++)
    {
      System.out.print(langNames[l] + "\t");
      for (int i =0; i < entries.size(); i++)
        System.out.print(entries.get(i).get(l) + "\t");
      System.out.println();
    }
  }
}
