package pty.smc.test;
import goblin.Taxon;

import java.io.*;
import java.util.*;


import nuts.util.Arbre;

import pty.RootedTree;
import pty.smc.PartialCoalescentState;

public class CountingDevice
{
  public static int countLeavesRelabelling(PartialCoalescentState p)
  {
    return countLeavesRelabelling(p.getFullCoalescentState());
  }
  
  public static int countLeavesRelabelling(RootedTree c)
  {
    throw new RuntimeException();
//    Arbre<Taxon> a = c.topology();
//    Map<Arbre<Taxon>,Set<Taxon>> descMap = Arbre.leavesMap(a);
//    System.out.println (descMap);
//    
//    int result = 1;
//    for (Arbre<Taxon> node : descMap.keySet())
//      if (!node.isLeaf())
//      {
//        if (node.getChildren().size() != 2)
//          throw new RuntimeException();
//        final Arbre<Taxon> 
//          l = node.getChildren().get(0),
//          r = node.getChildren().get(1);
//        final boolean bothChildrenLeaves = 
//          l.isLeaf() && r.isLeaf();
//        final int leftDesc = descMap.get(l).size(),
//                  rightDesc= descMap.get(r).size();
//                
//        if (!bothChildrenLeaves)
//          result *= Arithmetic.binomial(rightDesc+leftDesc, rightDesc);
//      }
//    return result;
  }
} 
