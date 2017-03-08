package pty.smc.test;
import java.io.*;
import java.util.*;

import pty.UnrootedTree;

import nuts.io.IO;

import ma.newick.NewickParser;

public class RemoveUselessTopNode
{

  /**
   * @param args
   */
  public static void main(String[] args)
  {
    UnrootedTree nct = UnrootedTree.fromNewick(new File("data/artif/rerootedInit.newick"));
    System.out.println(nct);
  }

}
