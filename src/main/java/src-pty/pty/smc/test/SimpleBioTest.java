package pty.smc.test;
import java.awt.geom.GeneralPath;
import java.io.*;
import java.util.*;

import nuts.lang.ArrayUtils;
import nuts.math.Fct;
import nuts.util.Arbre;
import nuts.util.Arbre.ArbreMap;

import pepper.Encodings;
import pty.eval.SymmetricDiff;
import pty.smc.PartialCoalescentState;
import pty.smc.ParticleFilter;
import pty.smc.ParticleKernel;
import pty.smc.PriorPostKernel;
import pty.smc.PriorPriorKernel;
import pty.smc.ParticleFilter.ParticleMapperProcessor;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;
import pty.smc.models.LikelihoodModelCalculator;

import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.Option;
import fig.exec.Execution;
import goblin.CognateId;
import goblin.Taxon;

import ma.BalibaseCorpus;
import ma.BioCorpus;
import ma.GeneratedCorpus;
import ma.MultiAlignment;
import ma.BalibaseCorpus.BalibaseCorpusOptions;
import ma.GeneratedCorpus.GeneratedCorpusOptions;
import ma.MultiAlignment.SequenceCoordinate;

public class SimpleBioTest implements Runnable
{
  @Option public Random rand = new Random(1);
  @Option public int N = 10;
  @Option public boolean useGenerated = true;
  @Option public double probabilityToMarkObservationAsMissing = 0.5;
  @Option public boolean usePriorPost = true;
  
  private BioCorpus corpus;
  private BalibaseCorpusOptions bco;
  private GeneratedCorpusOptions gco;
  public SimpleBioTest(BalibaseCorpusOptions bco, GeneratedCorpusOptions gco)
  {
    this.gco = gco;
    this.bco = bco;
  }
  public Set<Set<Taxon>> inferTree(CognateId id)
  {
    ParticleMapperProcessor<PartialCoalescentState, Set<Set<Taxon>>> mbr =
      SymmetricDiff.createCladeProcessor();
    ParticleKernel<PartialCoalescentState> kernel = (usePriorPost ?
        new PriorPostKernel(getInitState(id)) :
        new PriorPriorKernel(getInitState(id)));
    ParticleFilter.bootstrapFilter(kernel, mbr, N, rand);
//    LogInfo.logs(mbr.counter);
    for (Set<Set<Taxon>> cClades : mbr.getCounter())
      System.out.println(SymmetricDiff.clades2arbre(cClades).deepToLispString()
          + "\t" + mbr.getCounter().getCount(cClades));
    return mbr.centroid(SymmetricDiff.CLADE_SYMMETRIC_DIFFERENCE);
  }
  private PartialCoalescentState getInitState(CognateId id)
  {
    throw new RuntimeException();
//    MultiAlignment ma = corpus.getMultiAlignment(id);
//    Map<Language,List<Integer>> fullColumns = getFullColumns(corpus.getType().getEncodings(),ma);
//    List<Language> leafNames = new ArrayList<Language>();
//    List<LikelihoodModelCalculator> leaves = new ArrayList<LikelihoodModelCalculator>();
//    CTMC ctmc = CTMC.SimpleCTMC.fromSequenceType(fullColumns.values()
//        .iterator().next().size(),corpus.getType());
//    for (Language lang : fullColumns.keySet())
//    {
//      leafNames.add(lang);
//      leaves.add(DiscreteModelCalculator.observation(ctmc, toObservationArray(fullColumns.get(lang))));
//    }
//    return PartialCoalescentState.initialState(leaves,leafNames);
  }
  
  private int[] toObservationArray(List<Integer> list)
  {
    int [] result = new int[list.size()];
    for (int i = 0; i < result.length; i++)
      if (rand.nextDouble() < probabilityToMarkObservationAsMissing)
        result[i] = -1;
      else
        result[i] = list.get(i);
    return result;
  }

  public static void main(String [] args)
  {
    BalibaseCorpusOptions bco = new BalibaseCorpus.BalibaseCorpusOptions();
    GeneratedCorpusOptions gco = new GeneratedCorpusOptions();
    Execution.run(args, new SimpleBioTest(bco,gco), "bali", bco, "gen", gco);
  }
  
  public static Map<Taxon,List<Integer>> getFullColumns(Encodings enc, MultiAlignment ma)
  {
    Map<Taxon,List<Integer>> result = init(ma.getSequences().keySet());
    for (SequenceCoordinate sc : ma.eqClasses().representatives())
      if (!ma.isReference() || sc.isCoreBlock())
        if (ma.eqClasses().eqClass(sc).size() == ma.getSequences().size())
          process(ma.eqClasses().eqClass(sc), result, enc);
    return result;
  }

  private static void process(Set<SequenceCoordinate> eqClass,
      Map<Taxon, List<Integer>> result, Encodings enc)
  {
    for (SequenceCoordinate sc : eqClass)
      result.get(sc.getNodeIdentifier()).add(enc.char2PhoneId(sc.getCharValue()));
  }

  private static Map<Taxon, List<Integer>> init(Set<Taxon> keySet)
  {
    Map<Taxon, List<Integer>> result = new HashMap<Taxon, List<Integer>>();
    for (Taxon lang : keySet) result.put(lang,new ArrayList<Integer>());
    return result;
  }

  public void run()
  {
    if (useGenerated) corpus = new GeneratedCorpus(gco);
    else              corpus = new BalibaseCorpus (bco);
    for (CognateId id : corpus.intersectedIds())
      evaluate(id);
  }
  private void evaluate(CognateId id)
  {
    LogInfo.track("Evaluating " + id,true);
    Set<Set<Taxon>> recon = inferTree(id);
    Arbre<Taxon>  gold  = Arbre.tree2Arbre(corpus.getTopology(id)).preOrderMap(new ArbreMap<String,Taxon>() {
        @Override public Taxon map(Arbre<String> currentDomainNode) { return new Taxon(currentDomainNode.getContents()); }
      });
    Arbre<Taxon> reconTree = SymmetricDiff.clades2arbre(recon);
    Set<Set<Taxon>> goldClades = SymmetricDiff.clades(gold);
    LogInfo.logs("Gold=" + SymmetricDiff.clades2arbre(goldClades).deepToString());
    LogInfo.logs("Recon=" + reconTree.deepToString());
    LogInfo.logs("SymmetricCladeDiff(norm)="+SymmetricDiff.symmetricDifferenceSize(recon, goldClades)
        + "("+SymmetricDiff.normalizedSymmetricCladeDiff(reconTree, gold) +")");
    LogInfo.end_track();
  }
}
