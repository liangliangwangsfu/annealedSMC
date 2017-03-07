package pty.smc;
import goblin.Taxon;
import hmm.Param;
import hmm.ParamUtils;

import java.io.*;
import java.util.*;

import org.junit.Test;

import ev.io.PreprocessGutellData;
import ev.poi.processors.TreeDistancesProcessor;
import fig.basic.UnorderedPair;

import pepper.editmodel.OldMain;
import pty.io.Dataset;
import pty.io.Dataset.DatasetUtils;
import pty.smc.LazyParticleFilter.Eager2LazyAdaptor;
import pty.smc.LazyParticleFilter.LazyParticleKernel;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.ParticleFilter.ResamplingStrategy;
import pty.smc.models.CTMC;
import pty.smc.test.TestBrownianModel.KernelType;
import pty.smc.test.TestParticleNormalization.HMMParticleKernel;

import ma.MSAPoset;
import ma.SequenceType;
import nuts.math.TreeSumProd;
import nuts.math.TreeSumProd.HmmAdaptor;
import nuts.util.CollUtils.*;
import nuts.util.CollUtils;
import nuts.util.Counter;
import nuts.util.MathUtils;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public class LazyParticleFilterTest
{
  
  public static void main(String [] args)
  {
    new LazyParticleFilterTest().testTree();
  }
  
  @Test public void testTree()
  {
    ParticleFilter<LazyPCS> pc = new ParticleFilter<LazyPCS>();
    pc.N = 1000;
    pc.nThreads = 1;
    pc.resampleLastRound = true;
    pc.verbose = false;
    
    Random rand = new Random(1);
    MSAPoset msa =         
      PreprocessGutellData.randomDataSet(new File("/Users/bouchard/Documents/data/gutell/16S.3.alnfasta"), 1, 10, rand).get(0);
    
    
    Dataset dataset = DatasetUtils.fromAlignment(msa, SequenceType.RNA);
    
//    KernelType kernelType = KernelType.PRIOR_PRIOR;

    CTMC ctmc = CTMC.SimpleCTMC.dnaCTMC(dataset.nSites());
    PartialCoalescentState init = PartialCoalescentState.initFastState(false, dataset, ctmc);
    ParticleKernel<LazyPCS> ppk = new LazyPriorPrior(init);
    ParticleProcessor voidPro = new ParticleFilter.DoNothingProcessor();
    pc.sample(ppk,voidPro);
    
    
    System.out.println(pc.estimateNormalizer());
    Double newApprox = null;
    Counter<UnorderedPair<Taxon,Taxon>> ref = null;
    
    for (int i = 0; i < 4; i++)
    {
      LazyParticleKernel pk2 = new PriorPriorKernel(init);
//      Eager2LazyAdaptor adaptor = new Eager2LazyAdaptor(abject);
      LazyParticleFilter.ParticleFilterOptions options = new LazyParticleFilter.ParticleFilterOptions();
      options.nParticles = pc.N;
      options.rand = new Random(1);
      options.resampleLastRound = true;
      options.parallelizeFinalParticleProcessing = true;
      options.nThreads = 10 * i + 1;
      options.verbose = false;
      options.check();
      LazyParticleFilter lpf = new LazyParticleFilter(pk2, options);
      
      TreeDistancesProcessor tdp = new TreeDistancesProcessor();
      double zHat = lpf.sample(tdp);
      System.out.println("Approx2=" + zHat);
      
      if (ref == null)
        ref = tdp.getMeanDistances();
      
      Counter<UnorderedPair<Taxon, Taxon>> cur = tdp.getMeanDistances();
      for (UnorderedPair<Taxon, Taxon> p : CollUtils.union(ref.keySet(), cur.keySet()))
        MathUtils.checkClose(cur.getCount(p) ,ref.getCount(p));
//          throw new RuntimeException("" + cur.getCount(p) +" != " +  ref.getCount(p));
      
      if (newApprox == null)
        newApprox = zHat;
      if (newApprox != zHat)
        throw new RuntimeException();

    }
  }
  
  
  @Test public void testHMM()
  {
    Random rand = new Random(1);
    // generate a random set of params
    Param p = ParamUtils.randomUniParam(rand,10,10);
    System.out.println("Param:\n"+p);
    int [] obs = new int[20];
    System.out.println("Observation:" + Arrays.toString(obs));
    // compute normalization using particule filter
    HMMParticleKernel pk = new HMMParticleKernel(p, obs);
    ParticleProcessor voidPro = new ParticleFilter.DoNothingProcessor();
    ParticleFilter pf = new ParticleFilter();
    pf.N = 100;
    pf.resamplingStrategy = ResamplingStrategy.ALWAYS;
    pf.resampleLastRound = false;
    pf.sample(pk, voidPro);
    final double approxRef = pf.estimateNormalizer();
    System.out.println("Approx=" + approxRef);
    // exact
    ArrayList obsList = new ArrayList();
    for (int o : obs) obsList.add(o);
    HmmAdaptor adapt = new HmmAdaptor(p,obsList);
    TreeSumProd tsp = new TreeSumProd(adapt);
    final double exact = tsp.logZ();
    System.out.println("Exact="+exact);
    // new
    
    Double newApprox = null;
    for (int i = 0;i < 10; i++)
    {
      Eager2LazyAdaptor adaptor = new Eager2LazyAdaptor(pk);
      LazyParticleFilter.ParticleFilterOptions options = new LazyParticleFilter.ParticleFilterOptions();
      options.nParticles = pf.N;
//      options.minNDistinctParticles = pf.N / 10;
//      options.maxNParticles = pf.N * 100;
      options.resampleLastRound = false;
      options.nThreads = 10 * i + 1;
      options.verbose = true;
      options.check();
//        throw new RuntimeException();
      LazyParticleFilter lpf = new LazyParticleFilter(adaptor, options);
      double zHat = lpf.sample();
      System.out.println("Approx2=" + zHat);
      if (newApprox == null)
        newApprox = zHat;
      if (newApprox != zHat)
        throw new RuntimeException();
      
      MathUtils.threshold = 1.0;
      MathUtils.checkClose(newApprox, approxRef);
      MathUtils.checkClose(newApprox, exact);
      
      
    }
  }
}
