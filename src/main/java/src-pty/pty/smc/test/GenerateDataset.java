package pty.smc.test;
import java.io.*;
import java.util.*;

public class GenerateDataset
{
  public static void main(String [] args)
  {
    final Random rand = new Random(1);
    final double scale = 0.01;
    final double root = 0.5;
    final int nsites = 10000;
    double []
      isolate = new double[nsites],
      lang1   = new double[nsites],
      lang2   = new double[nsites],
      lang3   = new double[nsites],
      lang4   = new double[nsites],
      lang5 = new double[nsites],
      lang6 = new double[nsites];
      
    for (int i = 0; i < nsites; i++)
    {
      isolate[i] = fix(root + rand.nextGaussian()* Math.sqrt(scale * 5.0));
      double insideNode1 = root + rand.nextGaussian() * Math.sqrt(scale * 1.0);
      double insideNode2 = insideNode1 + rand.nextGaussian() * Math.sqrt(scale * 1.0);
      double insideNode3 = insideNode1 + rand.nextGaussian() * Math.sqrt(scale * 2.0);
      
      double insideNode4 = insideNode2 + rand.nextGaussian() * Math.sqrt(scale * 1.0);
      double insideNode5 = insideNode2 + rand.nextGaussian() * Math.sqrt(scale*1.5);
      
      lang1[i] = fix(insideNode4 + rand.nextGaussian() * Math.sqrt(scale * 1.0));
      lang2[i] = fix(insideNode4 + rand.nextGaussian() * Math.sqrt(scale * 1.0));
      lang3[i] = fix(insideNode5 + rand.nextGaussian() * Math.sqrt(scale * 0.5));
      lang4[i] = fix(insideNode5 + rand.nextGaussian() * Math.sqrt(scale * 0.5));
      
      lang5[i] = fix(insideNode3 + rand.nextGaussian() * Math.sqrt(scale * 1.0));
      lang6[i] = fix(insideNode3 + rand.nextGaussian() * Math.sqrt(scale * 1.0));
      
    }
    
    try {
    BufferedWriter bw = new BufferedWriter (new FileWriter("simulated.txt"));
    bw.write("isolate\t"); for (int i = 0; i < nsites; i++) bw.write(""+isolate[i]+"\t");
    bw.write("\nlang1\t");   for (int i = 0; i < nsites; i++) bw.write(""+lang1[i]+"\t");
    bw.write("\nlang2\t");   for (int i = 0; i < nsites; i++) bw.write(""+lang2[i]+"\t");
    bw.write("\nlang3\t");   for (int i = 0; i < nsites; i++) bw.write(""+lang3[i]+"\t");
    bw.write("\nlang4\t");   for (int i = 0; i < nsites; i++) bw.write(""+lang4[i]+"\t");
    bw.write("\nlang5\t");   for (int i = 0; i < nsites; i++) bw.write(""+lang5[i]+"\t");
    bw.write("\nlang6\t");   for (int i = 0; i < nsites; i++) bw.write(""+lang6[i]+"\t");
    bw.close ();
    } catch (Exception e) {
    	e.printStackTrace ();
    }
  }

  private static double fix(double d)
  {
    if (d < 0)
    {
      System.err.println("Less than zero");
      return 0;
    }
    else if (d > 1)
    {
      System.err.println("More than one");
      return 1;
    }
    else return d;
  }
  
}
