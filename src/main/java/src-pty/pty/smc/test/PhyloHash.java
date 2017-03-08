package pty.smc.test;
import java.io.*;
import java.util.*;

import nuts.util.CollUtils;

import pty.RootedTree;
import pty.smc.PartialCoalescentState;

public class PhyloHash
{

  public static boolean equal(PartialCoalescentState p1, PartialCoalescentState p2)
  {
    return p1.getPhylo().equals(p2.getPhylo());
  }
  
  public static int hashCode(PartialCoalescentState p)
  { 
    return p.getPhylo().hashCode();
  }
  
  public static class Phylo<T>
  {
    private final T field;
    private final Set<Phylo<T>> children;
    public Phylo(T field, Set<Phylo<T>> children)
    {
      this.field = field;
      this.children = CollUtils.archive(children);
    }
    @Override
    public int hashCode()
    {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((children == null) ? 0 : children.hashCode());
      result = prime * result + ((field == null) ? 0 : field.hashCode());
      return result;
    }
    @Override
    public boolean equals(Object obj)
    {
      if (this == obj)
        return true;
      if (obj == null)
        return false;
      if (getClass() != obj.getClass())
        return false;
      Phylo other = (Phylo) obj;
      if (children == null)
      {
        if (other.children != null)
          return false;
      } else if (!children.equals(other.children))
        return false;
      if (field == null)
      {
        if (other.field != null)
          return false;
      } else if (!field.equals(other.field))
        return false;
      return true;
    }
  }
}