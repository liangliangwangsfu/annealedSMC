package pty.io;
import java.io.*;
import java.util.*;

import pty.UnrootedTree;

public class ShowTree
{
  public static void main(String [] args)
  {
    UnrootedTree nct = UnrootedTree.fromNewick(new File("data/world-soft-hgdp.txt"));
    System.out.println(nct);
  }
}
