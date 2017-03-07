package pty.smc;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import pty.io.Dataset;
import pty.io.HGDPDataset;
import pty.io.WalsDataset;
import pty.io.Dataset.DatasetType;
import pty.smc.models.BrownianModel;
import pty.smc.models.BrownianModelCalculator;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;
import pty.smc.test.PhylipParser;
import pty.smc.test.TestJointModel;

import fig.basic.IOUtils;
import fig.exec.Execution;
import goblin.Taxon;
import ma.newick.NewickParser;
import ma.newick.ParseException;
import nuts.util.Arbre;
import nuts.util.Tree;

public class TestLikelihood implements Runnable {

	Arbre<String> inputArbre;
	
	
	private void init (Dataset data) {
		BrownianModel bm = new BrownianModel(1, 1);
		List<Taxon> leafNames = new ArrayList<Taxon>();
		List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
		Map<Taxon, double[][]> observations = data.observations();
		int n = observations.keySet().size();
		for (Taxon lang : observations.keySet())
		{
			
		      leafNames.add(lang);
		      double [][] cObs = observations.get(lang);
		      double [] converted = new double[cObs.length];
		      for (int i = 0; i < converted.length;i++)
		        converted[i] = cObs[i][0];
		      
		      
		      leaves.add( BrownianModelCalculator.observation (converted, bm, false));
		}
		
		while (n>1) {
			List<LikelihoodModelCalculator> tmp = new ArrayList<LikelihoodModelCalculator> ();
			
			for (int i = 0 ; i < n; i=i+2) {
				tmp.add(leaves.get(i).combine (leaves.get(i), leaves.get(i+1),1,1, false));
			}
			leaves = tmp;
			n = leaves.size ();
		}
		
		System.out.println ("Log l ="+leaves.get(0).logLikelihood());
		double llnew = - (Math.log (6) + 0.5*Math.log(7) + 1.5 *Math.log (Math.PI));
		double llold = - (Math.log (96) + 0.5*Math.log(5) + 3.5 *Math.log (Math.PI));
		double lz1 = -0.5 * Math.log (4*Math.PI);
		double lz2 = -0.5 * Math.log (6*Math.PI);
		double lz3 = -0.5 * Math.log(7*Math.PI);
		double ll1 = lz1;
		double ll2 = lz2 + 2*lz1;
		double ll3 = lz3 + 2*lz2;
		
		System.out.println (ll1 + "\t" + ll2 + "\t" + ll3);
		System.out.println ("Log l  =" + llnew + "\t" + llold);  
	}
	
	
//use TestJointModel.initGeneState()
//	private PartialCoalescentState initGeneState (String filename) { 
//		PhylipParser pp = new PhylipParser (filename);
//		int nleaves =  pp.getNumberOfLeaves();
//		int nsites = pp.getNumberofSites();
//		double[][] observations = pp.getDoubleObservations();
//		ArrayList<Language> leafNames = pp.getLeafNames();
//		
//		double rateScalar = 0.5;
//		double[][] rate = {{-rateScalar,+rateScalar}, {+rateScalar,-rateScalar}} ;
//		CTMC ctmc = new CTMC.SimpleCTMC (rate, nsites);
//		ArrayList<LikelihoodModelCalculator> dmcList = new ArrayList <LikelihoodModelCalculator> ();
//		for (int i = 0 ; i < nleaves; i++) {
//			dmcList.add( DiscreteModelCalculator.observation (ctmc, observations[i]));
//		}
//		PartialCoalescentState pcs = PartialCoalescentState.initialState( dmcList, leafNames);
//		return pcs;		
//		
//	}

	
	public void computeLogLikelihood () {
		
	}
	
  public static void main (String args []){
    Execution.run(args, new TestLikelihood(),
        "hgdp", HGDPDataset.class);
  }

  public void run()
  {
	  Dataset data = DatasetType.HGDP.loadDataset();
	  init(data);
  }
}
