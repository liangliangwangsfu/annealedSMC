package pty.eval;
import java.io.*;
import java.util.*;

import nuts.util.CollUtils;

public class CladeOverlap
{
  public static <T> double cladeOverlap(Set<Set<T>> clades1, Set<Set<T>> clades2)
  {
    double sum = 0.0;
    for (Set<T> clade1 : clades1)
      for (Set<T> clade2 : clades2)
        sum += _cladeOverlap(clade1, clade2);
    return sum;
  }
  public static <T> double _cladeOverlap(Set<T> clades1, Set<T> clades2)
  {
    final double interSize = CollUtils.inter(clades1,clades2).size();
    return interSize / (clades1.size() + clades2.size() - interSize);
//    return 1.0 / ( (clades1.size() + clades2.size()) / CollUtils.inter(clades1,clades2).size() - 1);
  }
}
