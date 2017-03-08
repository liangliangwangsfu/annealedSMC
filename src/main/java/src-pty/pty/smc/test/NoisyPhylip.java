package pty.smc.test;
import java.io.*;
import java.util.*;

import pty.io.WalsAnn;

public class NoisyPhylip
{
  public static void main(String [] args)
  {
    Random rand = new Random(1);
    int nSites = 10000;
    int nPops = 14;
    System.out.println(" " + nPops + " " + nSites);
    for (int p = 0; p < nPops; p++)
    {
      System.out.print(WalsAnn.cleanForPhylip("Pop"+p));
      for (int s =0 ; s < nSites;s++)
        System.out.print(rand.nextInt(2)+"\t");
      System.out.println();
    }
  }
}
