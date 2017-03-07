package pty.io;
import fig.basic.NumUtils;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import nuts.lang.ArrayUtils;
import nuts.math.GMFct;
import nuts.math.TreeSumProd;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import pty.RootedTree;
import pty.Observations;
import pty.learn.DiscreteBP;
import pty.smc.PartialCoalescentState;
import pty.smc.models.CTMC;

public class LeaveOneOut
{
  public static double loo(PartialCoalescentState pcs)
  {
    return loo(pcs.getFullCoalescentState(), pcs.getCTMC(), pcs.getObservations());
  }
  public static double loo(
      RootedTree tree,
      CTMC ctmc,
      Observations observations)
  {
    SummaryStatistics stats = new SummaryStatistics();
    // loop over modern species
    for (Taxon lang : observations.observations().keySet())
      for (int site = 0; site < observations.nSites(); site++)
        if (isKnown(observations,lang,site))
        {
          final GMFct<Taxon> post = DiscreteBP.posteriorMarginalTransitions(
              tree, ctmc, observations, site, lang);
          int prediction = -1;
          double max = Double.NEGATIVE_INFINITY, current;
          for (int ch = 0; ch < observations.nCharacter(site); ch++)
            if ((current = post.get(lang, ch)) > max)
            {
              prediction = ch;
              max = current;
            }
          final int truth = getObs(observations,lang,site);
          stats.addValue((prediction == truth ? +1 : 0));
        }
    return stats.getMean();
  }
  // hacks :(
  private static int getObs(Observations observations, Taxon lang, int site)
  {
    int result = -1;
    final double [] array = observations.observations().get(lang)[site];
    for (int i = 0; i < array.length; i++)
    {
      final double cur = array[i]; 
      if (cur == 1.0)
      {
        if (result == -1) result = i;
        else throw new RuntimeException();
      }
    }
    return result;
  }
  private static boolean isKnown(Observations observations, Taxon lang, int site)
  {
    double sum = ArrayUtils.sum(observations.observations().get(lang)[site]);
    if (sum == 1.0) return true;
    else            return false;
  }
}
