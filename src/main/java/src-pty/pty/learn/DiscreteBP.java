package pty.learn;
import fig.basic.LogInfo;
import fig.basic.Pair;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import pty.RootedTree;
import pty.Observations;
import pty.eval.SymmetricDiff;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.PriorPriorKernel;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.ParticleFilter.ParticleProcessor;
import pty.smc.models.CTMC;
import pty.smc.test.SymmetryTest;
import pty.smc.test.SymmetryTest.MeanHeight;

import nuts.math.GMFct;
import nuts.math.GMFctUtils;
import nuts.math.Graph;
import nuts.math.HashGraph;
import nuts.math.TabularGMFct;
import nuts.math.TreeSumProd;
import nuts.util.Arbre;
import nuts.util.Tree;

/**
 * Adaptor for doing discrete Belief Prop. on tree structure
 * @author bouchard
 *
 */
public class DiscreteBP
{
//  public static GMFct<Language> posteriorMarginalTransitions(PartialCoalescentState pcs, int site)
//  {
////    CTMC ctmc = pcs.getCTMC();
////    Pair<Arbre<Language>,Map<Language,Double>> p = pcs.getArbreAndBranchLengths();
////    Arbre<Language> a = p.getFirst(); Map<Language,Double> bls = p.getSecond(); 
////    Map<Language,double[]> obs = pcs.getObservations(site);
//    return posteriorMarginalTransitions(a,bls,ctmc,obs,site);
//  }
  
  
  public static GMFct<Taxon> posteriorMarginalTransitions(
      PartialCoalescentState pcs,
      final int site)
  {
    return TreeSumProd.computeMoments(toGraphicalModel(pcs.getFullCoalescentState(),
        pcs.getCTMC(),
        pcs.getObservations(),site));
  }
  
  public static double dataLogLikelihood(
      RootedTree state,
      CTMC ctmc,
      Observations observations,
      final int site)
  {
    return new TreeSumProd<Taxon>(toGraphicalModel(state,ctmc,observations,site)).logZ();
  }
  
  public static double dataLogLikelihood(
      RootedTree state,
      CTMC ctmc,
      Observations observations)
  {
    double sum = 0.0; 
    for (int s = 0; s < ctmc.nSites(); s++)
      sum += dataLogLikelihood(state,ctmc,observations, s);
    return sum;
  }
  
  public static GMFct<Taxon> posteriorMarginalTransitions(
      RootedTree state,
      CTMC ctmc,
      Observations observations,
      final int site)
  { return posteriorMarginalTransitions(state, ctmc, observations, site, null); }
  
  public static GMFct<Taxon> posteriorMarginalTransitions(
      RootedTree state,
      CTMC ctmc,
      Observations observations,
      final int site,
      Taxon languageToHeldout)
  {
    return TreeSumProd.computeMoments(toGraphicalModel(state,ctmc,observations,site, languageToHeldout));
  }
  
  public static GMFct<Taxon> toGraphicalModel(
      RootedTree state,
      CTMC ctmc,
      Observations observations,
      final int site)
  { return toGraphicalModel(state,ctmc,observations,site,null); }
  public static GMFct<Taxon> toGraphicalModel(
      RootedTree state,
      CTMC ctmc,
      Observations observations,
      final int site,
      Taxon languageToHeldout)
  {
    final int nCharacters = ctmc.nCharacter(site);
    // create the graph
    Tree<Taxon> t = Arbre.arbre2Tree(state.topology());
    Graph<Taxon> g = new HashGraph<Taxon>(t);
    // create the transition structure
    Map<Taxon,Integer> rvRange = new HashMap<Taxon,Integer>();
    for (Taxon l : state.topology().nodeContents()) rvRange.put(l,nCharacters);
    TabularGMFct<Taxon> pot = GMFctUtils.ones(new TabularGMFct<Taxon>(g, rvRange));
    // handle transitions
    for (Arbre<Taxon> node : state.topology().nodes())
      if (!node.isRoot())
      {
        Taxon cur = node.getContents(), par = node.getParent().getContents();
        double [][] prs = ctmc.getTransitionPr(site, state.branchLengths().get(node.getContents()));
//        LogInfo.logsForce("tgm prs:" + Arrays.deepToString(prs));
        for (int ts = 0; ts < nCharacters; ts++)
          for (int bs = 0; bs < nCharacters; bs++)
            pot.set(par, cur, ts, bs, prs[ts][bs]);
      }
    // handle initial distribution
    Taxon root = state.topology().getContents();
    double [] initDist = ctmc.getInitialDistribution(site);
//    LogInfo.logsForce("tgm initD:" + Arrays.toString(initDist));
    for (int c = 0; c < nCharacters; c++)
      pot.set(root, c, initDist[c]);
    // handle observations
    Map<Taxon,double[][]> obs = observations.observations();
    for (Taxon lang : obs.keySet())
      for (int c = 0; c < nCharacters; c++)
      {
//        LogInfo.logsForce("tgm obs:lang=" + lang + ",c=" + c + ",value=" + obs.get(lang)[c]);
        if (lang.equals(languageToHeldout))
          pot.set(lang, c, 1.0);
        else
          pot.set(lang, c, obs.get(lang)[site][c]);
      }
    return pot;
  }
//  public static GMFct<Language> toGraphicalModel(PartialCoalescentState pcs, int site)
//  {
//    CTMC ctmc = pcs.getCTMC();
//    Pair<Arbre<Language>,Map<Language,Double>> p = pcs.getArbreAndBranchLengths();
//    Arbre<Language> a = p.getFirst(); Map<Language,Double> bls = p.getSecond(); 
//    Map<Language,double[]> obs = pcs.getObservations(site);
//    return toGraphicalModel(a,bls,ctmc,obs,site);
//    CTMC ctmc = pcs.getCTMC();
//    final int nCharacters = ctmc.nCharacter(site);
//    Pair<Arbre<Language>,Map<Language,Double>> p = pcs.getArbreAndBranchLengths();
//    Arbre<Language> a = p.getFirst(); Map<Language,Double> bls = p.getSecond();
//    // create the graph
//    Tree<Language> t = Arbre.arbre2Tree(a);
//    Graph<Language> g = new Graph.HashGraph<Language>(t);
//    // create the transition structure
//    Map<Language,Integer> rvRange = new HashMap<Language,Integer>();
//    for (Language l : a.nodeContents()) rvRange.put(l,nCharacters);
//    TabularGMFct<Language> pot = GMFctUtils.ones(new TabularGMFct<Language>(g, rvRange));
//    // handle transitions
//    for (Arbre<Language> node : a.nodes())
//      if (!node.isRoot())
//      {
//        Language cur = node.getContents(), par = node.getParent().getContents();
//        double [][] prs = ctmc.getTransitionPr(site, bls.get(node.getContents()));
////        LogInfo.logsForce("tgm prs:" + Arrays.deepToString(prs));
//        for (int ts = 0; ts < nCharacters; ts++)
//          for (int bs = 0; bs < nCharacters; bs++)
//            pot.set(par, cur, ts, bs, prs[ts][bs]);
//      }
//    // handle initial distribution
//    Language root = a.getContents();
//    double [] initDist = ctmc.getInitialDistribution(site);
////    LogInfo.logsForce("tgm initD:" + Arrays.toString(initDist));
//    for (int c = 0; c < nCharacters; c++)
//      pot.set(root, c, initDist[c]);
//    // handle observations
//    Map<Language,double[]> obs = pcs.getObservations(site);
//    for (Language lang : obs.keySet())
//      for (int c = 0; c < nCharacters; c++)
//      {
////        LogInfo.logsForce("tgm obs:lang=" + lang + ",c=" + c + ",value=" + obs.get(lang)[c]);
//        pot.set(lang, c, obs.get(lang)[c]);
//      }
//    return pot;
//  }
  public static void main(String [] args)
  {
//    Random rand = new Random(1);
//    int N = 10;
//    SaveParticles<PartialCoalescentState> sp = new SaveParticles<PartialCoalescentState>();
//    PartialCoalescentState pcs = SymmetryTest.getInitState(4,1);
//    PriorPriorKernel kernel = new PriorPriorKernel(pcs);
//    ParticleFilter.bootstrapFilter(kernel, sp, N, rand);
//    PartialCoalescentState fcs = sp.states().iterator().next();
//    System.out.println(fcs.getArbreAndBranchLengths().getFirst().deepToString());
//    System.out.println(GMFctUtils.toString(posteriorMarginalTransitions(fcs,0)));
  }
}
