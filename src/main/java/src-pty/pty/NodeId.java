package pty;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;


public class NodeId// implements Serializable, Comparable<NodeId>
{
//  private static Map<String,NodeId> canonicalMap = new HashMap<String,NodeId>();
//  public static NodeId getInstance(String langStr)
//  {
//    NodeId result = canonicalMap.get(langStr);
//    if (result == null)
//    {
//      result = new NodeId(langStr);
//      result.isCan = true;
//      synchronized (canonicalMap) { canonicalMap.put(langStr, result); }
//    }
//    return result;
//  }
//  
//  private static final long serialVersionUID = 1L;
//  private final String string;
//  private boolean isCan = false;
//  private NodeId(String string) 
//  { 
//    if (string == null)
//      throw new RuntimeException("Name of lang should be nontrivial");
//    this.string = string; 
//  }
//  @Override
//  public boolean equals(Object obj)
//  {
//    if (!isCan) throw new RuntimeException();
//    return super.equals(obj);
//  }
//  @Override
//  public int hashCode()
//  {
//    if (!isCan) throw new RuntimeException();
//    return super.hashCode();
//  }
//  @Override
//  public String toString()
//  {
//    return string;
//  }
//  public int compareTo(NodeId arg0)
//  {
//    return this.string.compareTo(arg0.string);
//  }
  
}
