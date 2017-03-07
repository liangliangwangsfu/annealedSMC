package pty.smc.test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import ma.newick.NewickParser;
import ma.newick.ParseException;
import nuts.maxent.SloppyMath;
import nuts.util.Arbre;
import nuts.util.Counter;
import pty.RootedTree;
import pty.Train;
import pty.eval.Purity;
import pty.eval.SymmetricDiff;
import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.learn.CTMCLoader;
import pty.smc.ConditionalPriorPriorKernel;
import pty.smc.ConstrainedKernel;
import pty.smc.MapLeaves;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PostPostKernel;
import pty.smc.PostPostKernelA;
import pty.smc.PriorPriorKernel;
import pty.smc.PriorPostKernel;
import pty.smc.ParticleFilter.ForkedProcessor;
import pty.smc.ParticleFilter.PCSHash;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.models.BrownianModel;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import fig.basic.IOUtils;
import fig.basic.LogInfo;
import fig.basic.Option;
import fig.basic.Pair;
import fig.exec.Execution;
import goblin.BayesRiskMinimizer;
import goblin.Taxon;


public class TestBrownianModel implements Runnable {
//	@Option public String genefile;  // Use HGDP.path instead
//	@Option public String mapfile;
//	@Option public int gibbsIterations = 5;
	@Option public static final boolean test = false;
	
  public static ParticleFilter<PartialCoalescentState> pf = new ParticleFilter<PartialCoalescentState>();
	public static final Random rand= new Random(1);
	
	@Option public KernelType kernelType = KernelType.POST_POST;
	@Option public String constraintsFile = "";
	@Option public boolean usePriorPost = false; // used only in constrained kernel
	@Option public double variance = 1;
	@Option public boolean linesearch = true;
	
	public static enum KernelType 
	{ 
	  PRIOR_POST {
      @Override public ParticleKernel<PartialCoalescentState> load(
          PartialCoalescentState initial,
          String constraintFile)
      {
        return new PriorPostKernel(initial);
      }
    }, 
    POST_POST {
      @Override public ParticleKernel<PartialCoalescentState> load(
          PartialCoalescentState initial,
          String constraintFile)
      {
          return new PostPostKernel(initial);
//        return new PostPostKernelA(initial);
      }
    }, 
	  PRIOR_PRIOR {
      @Override public ParticleKernel<PartialCoalescentState> load(
          PartialCoalescentState initial, 
          String constraintFile)
      {
        return new PriorPriorKernel(initial);
      }
    }, 
    PRIOR_POST2 {
      @Override public ParticleKernel<PartialCoalescentState> load(
          PartialCoalescentState initial,
          String constraintFile)
      {
        ConstrainedKernel.usesPriorPost = true;
        return new ConstrainedKernel(initial, new HashSet<Set<Taxon>>());
      }
    }, 
	  CONSTRAINED {
      @Override public ParticleKernel<PartialCoalescentState> load(
          PartialCoalescentState initial, String constraintFile)
      {
        try
        {
          final NewickParser np = new NewickParser(IOUtils.openInHard(constraintFile));
          final Arbre<Taxon> a = Taxon.LanguageUtils.convert(Arbre.tree2Arbre(np.parse()));
          return new ConstrainedKernel(initial, SymmetricDiff.clades(a)); //FromUnrooted(a));
        } 
        catch (ParseException e) { throw new RuntimeException(e); }
      }
    };
	  public abstract ParticleKernel<PartialCoalescentState> load(PartialCoalescentState initial, String constraintFile);
	} 
	
	PartialCoalescentState genepcs, langpcs;
	public static void main (String args []){
    Execution.monitor = true;
    Execution.makeThunk = false;
    Execution.create = true;
    Execution.useStandardExecPoolDirStrategy = true;
		Execution.run(args, new TestBrownianModel(),
        "hgdp", HGDPDataset.class,
        "kernel", ConstrainedKernel.class,
        "filter", pf,
        "pcs", PartialCoalescentState.class,
        "ppk", PriorPriorKernel.class);
	}
	

	public static PartialCoalescentState initGeneState (double variance) {	  
	  Dataset data = DatasetType.HGDP.loadDataset();
	  BrownianModel bm = new BrownianModel(data.nSites(), variance);
	  return PartialCoalescentState.initState(data, bm,false);

	}
	
	
	@SuppressWarnings("unchecked")
  public void sample (PartialCoalescentState initGeneState) {
	//	ParticleKernel<PartialCoalescentState> ppk =
		//	new PriorPriorKernel (initGeneState);
		
		ParticleKernel<PartialCoalescentState> ppk = kernelType.load(initGeneState, constraintsFile);
		
		 ParticleFilter.ParticleMapperProcessor<PartialCoalescentState,PartialCoalescentState> mapd = 
		      ParticleFilter.ParticleMapperProcessor.saveCoalescentParticlesProcessor();
		 ParticleFilter.ParticleMapperProcessor<PartialCoalescentState,Set<Set<Taxon>>> mbr =
       SymmetricDiff.createCladeProcessor();
		 PCSHash hashProcessor = new PCSHash();
		 ForkedProcessor<PartialCoalescentState> forked = new ForkedProcessor<PartialCoalescentState>(mapd, mbr, hashProcessor);
		 pf.sample ( ppk , forked);
		 PartialCoalescentState pcs = mapd.map();
		 
		 if (test) {
			 ArrayList<PartialCoalescentState> samples = (ArrayList<PartialCoalescentState>)pf.getSamples();
			 double[] logWeights = pf.getLogWeights();
			 HashMap <Integer, Double> pcs2weightmap = new HashMap <Integer, Double>();
			 HashMap <Integer, Integer> pcs2countmap = new HashMap <Integer, Integer> ();
			 
			 for (int i = 0 ; i < samples.size(); i++) {
				 PartialCoalescentState p = samples.get(i);
				 double w = p.logLikelihood();
				 //System.out.println (p );
				 //System.out.println (p.hashCode());
				 pcs2weightmap.put(PhyloHash.hashCode(p), w);
				 pcs2countmap.put(PhyloHash.hashCode(p), CountingDevice.countLeavesRelabelling(p));
			 }
			 
			 double[] tmpw = new double[pcs2weightmap.keySet().size()];
			 int j = 0 ;
			 int numTrees =0 ;
			 for (Integer i:pcs2weightmap.keySet()) {
				 tmpw[j] = pcs2weightmap.get(i) + Math.log(pcs2countmap.get(i));
				 numTrees += pcs2countmap.get(i);
				 j++;
			 }
			 double marginallogl = SloppyMath.logAdd(tmpw) - Math.log(numTrees);
			 System.out.println ("******Test*****");
			 System.out.println ("Unique topologies = " + pcs2weightmap.keySet().size());			 			 
			 System.out.println ("Marginal = " + marginallogl + "?=" + pf.estimateNormalizer());
			 System.out.println ("Total num trees = " + numTrees);
			 
			 j = pcs.getLeaves();
			 int count = 1;
			 while ( j > 2) {
				count *= (j*(j-1)/2);
				j--;
			 }
			 System.out.println ("Total num trees (analytical) = " + count);
			 
		 }
		 
		 RootedTree map = pcs.getFullCoalescentState();
		 String prefix = "map-tree.variance-" + variance;
		 Train.outputTree (map.topology(), prefix, map.branchLengths());
		 Train.outputTree ( SymmetricDiff.clades2arbre(mbr.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE)), "mbr-tree");
		 LogInfo.logs("Log likelihood of map-tree:" + pcs.logLikelihood() + " with variance = " + variance);
		 LogInfo.logs("Marginal log likelihood:" + pf.estimateNormalizer() + " with variance = " + variance);
		 LogInfo.logs("Height of map tree: " + pcs.topHeight());
		 LogInfo.logs("Hash=" + hashProcessor.getHash());

	}

	public static final double logNormalDensity (final double x, final double mean, final double var) {
	    return -0.5*(x-mean)*(x-mean)/var -0.5*Math.log(2*Math.PI *var);
	}
	  
	public void run () {
		if (linesearch) {
			variance = 1e5;
			for (int i = 0 ; i<=10 ; i++, variance/=10) {
				PartialCoalescentState geneState = initGeneState (variance);
				sample (geneState);
			}
		} else { 
			PartialCoalescentState geneState = initGeneState (variance);
			sample (geneState);
		}
		
		double lz1 = logNormalDensity(0, 0, 8.5);
		double lz2 = logNormalDensity(0, 0, 3 + 72/17);
		double ll = lz2 + 2 * lz1;
		
//		double lz1 = logNormalDensity(0, 0, 2);
//	double lz2 = logNormalDensity(0, 0, 3.5);
//		double ll = lz2 + lz1;
		
		System.out.println ("ll = " + ll + "\tlz1 = " + lz1);
		
		double lz = logNormalDensity( 0, 0, 2);
		System.out.println ("Height 1 = " + lz);
		lz = logNormalDensity( 0, 0, 3);
		System.out.println ("Height 2 = " + lz);
		lz = logNormalDensity( 0, 0, 4);
		System.out.println ("Height 3 = " + lz);
		lz = logNormalDensity( 0, 0, 5);
		System.out.println ("Height 1 = " + lz);
		
	}	
	
}
