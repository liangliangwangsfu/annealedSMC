package pty.smc.test;
import goblin.CognateId;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.w3c.dom.ranges.RangeException;

import pty.eval.SymmetricDiff;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.PriorPriorKernel;
import pty.smc.PriorPostKernel;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import pty.smc.models.NoLikelihoodModel;

public class SymmetryTest
{
  public static final int N = 5;
  public static final Random rand= new Random(1);
  public static void main(String [] args)
  {
	/* 
	 * 
	 */
    ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> mbr =
      SymmetricDiff.createCladeProcessor();
    PartialCoalescentState pcs = getInitState(10,10);
//    PriorPriorKernel kernel = new PriorPriorKernel(pcs);
    PriorPostKernel kernel = new PriorPostKernel(pcs);
    MeanHeight mh = new MeanHeight();
    ParticleFilter.bootstrapFilter(kernel, mbr, 
//      mh,
        
//        new ParticleFilter.ParticleProcessor<PartialCoalescentState>() {
//          public void process(PartialCoalescentState state, double weight)
//          {
//            System.out.println("Weight=" + weight + ",Sample=\n" + state.toString());
//          }
//        },
        N, rand);
    for (Set<Set<Taxon>> cClades : mbr.getCounter())
      System.out.println(SymmetricDiff.clades2arbre(cClades).deepToLispString()
          + "\t" + mbr.getCounter().getCount(cClades));
//    System.out.println("Mean h=" + mh.ss.getSum());
    
  }
  
  public static class MeanHeight implements ParticleProcessor<PartialCoalescentState>
  {
    public final SummaryStatistics ss = new SummaryStatistics(), temp = new SummaryStatistics();
    public void process(PartialCoalescentState state, double weight)
    {
      ss.addValue(weight*state.topHeight());
      temp.addValue(weight);
    }
    
  }
  
  public static PartialCoalescentState getInitState(int nLeaves, int nSites)
  {
    throw new RuntimeException();
//    List<Language> leafNames = new ArrayList<Language>();
//    List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
//    double [][] rate = {{-0.5,+0.5},{+0.5,-0.5}};
//    CTMC ctmc = new CTMC.SimpleCTMC(rate,nSites);
//    int [] states0 = new int[nSites];
//    int [] states1 = new int[nSites];
//    for (int s = 0; s < nSites; s++)
//      states1[s] = 1;
//    for (int i = 0; i < nLeaves; i++)
//    {
//      leafNames.add(new Language("" + i));
//      leaves.add(DiscreteModelCalculator.observation(ctmc, (i<2 ? states0 : states1)));
////      leaves.add(new NoLikelihoodModel());
//    }
//    return PartialCoalescentState.initialState(leaves,leafNames);
  }

}
