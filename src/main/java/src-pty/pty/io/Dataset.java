package pty.io;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import fig.basic.StrUtils;
import goblin.Taxon;

import java.io.*;
import java.util.*;

import ma.MSAPoset;
import ma.SequenceType;
import ma.MSAPoset.Column;
import nuts.util.Arbre;
import nuts.util.CollUtils;
import nuts.util.Indexer;

import pty.ObservationDimensions;
import pty.Observations;
import pty.smc.PartialCoalescentState.CoalescentNode;
import pty.smc.models.CTMC;
import pty.smc.models.DiscreteModelCalculator;


public interface Dataset extends Observations
{
  public static final int unknownCode = -1;
  
  public boolean hasReferenceClusters();
  public Map<Taxon,String> getReferenceClusters();
  
  public static enum DatasetType
  {
    WALS { @Override public Dataset loadDataset() { return WalsDataset.getPreprocessedCorpus(); }},
    HGDP { @Override public Dataset loadDataset() { return new HGDPDataset(); }},
    GEN  { @Override public Dataset loadDataset() { return new GeneratedDataset(); }};
    public abstract Dataset loadDataset();
  }
  
  public static class DatasetUtils
  {
    public static Dataset fromAlignment(File alignmentFile,SequenceType sequenceType)
    {
      // read alignment
      MSAPoset align = MSAPoset.parseAlnOrMsfFormats(alignmentFile);
      return fromAlignment(align, sequenceType);
    }
    
    public static Dataset fromAlignment(MSAPoset align,SequenceType sequenceType)
    {
      // supports missing data
      // WARNING: copied from DiscreteModelCalculator...
      // create the appropriate  encodings
      Indexer<Character> indexer = sequenceType.getEncodings().nonGapCharactersIndexer();

      
      
      Map<Taxon,int[]> data = new HashMap<Taxon,int[]>();
      for (Taxon t : align.sequences().keySet())
        data.put(t, new int[align.columns().size()]);
      int i = 0;
      for (Column c : align.linearizedColumns())
      {
        for (Taxon t : align.sequences().keySet())
          if (c.getPoints().keySet().contains(t))
          {
            Character lastChar = null;
            try
            {
              lastChar = align.charAt(c, t);
              int index = indexer.o2i(lastChar);
              data.get(t)[i] = index;
            }
            catch (Exception e)
            {
//              LogInfo.warning("Probably an unknown character:" + lastChar + "... replacing by UNK");
              data.get(t)[i] = Dataset.unknownCode;
            }
          }
          else
            data.get(t)[i] = Dataset.unknownCode;
        i++;
      }
//        else
//          throw new RuntimeException("For now, only gap-free alignments accepted");
      final int nSites = align.columns().size();
      
      
      final int nChars = indexer.size();
      final ObservationDimensions od = new ObservationDimensions() {
        @Override
        public int nSites()
        {
         return nSites;
        }
        @Override
        public int nCharacter(int site)
        {
          return nChars;
        }
      };
      final Map<Taxon,double[][]> converted = convert(data, od);
      return new Dataset() {
        @Override
        public Map<Taxon, String> getReferenceClusters() { throw new RuntimeException(); }
        @Override
        public boolean hasReferenceClusters()
        {
          return false;
        }
        @Override
        public Map<Taxon, double[][]> observations()
        {
          return converted;
        }
        @Override
        public int nCharacter(int site)
        {
          return nChars;
        }
        @Override
        public int nSites()
        {
          return nSites;
        }
      };
    }
    
    
    public static Dataset fromAlignment(MSAPoset align,SequenceType sequenceType, int repeats)
    {
      // supports missing data
      // WARNING: copied from DiscreteModelCalculator...
      // create the appropriate  encodings
      Indexer<Character> indexer = sequenceType.getEncodings().nonGapCharactersIndexer();
      
      Map<Taxon,int[]> data = new HashMap<Taxon,int[]>();
      for (Taxon t : align.sequences().keySet())
        data.put(t, new int[align.columns().size()]);
      int i = 0;
      for (Column c : align.linearizedColumns())
      {
        for (Taxon t : align.sequences().keySet())
          if (c.getPoints().keySet().contains(t))
          {
            Character lastChar = null;
            try
            {
              lastChar = align.charAt(c, t);
              int index = indexer.o2i(lastChar);
              data.get(t)[i] = index;
            }
            catch (Exception e)
            {
//              LogInfo.warning("Probably an unknown character:" + lastChar + "... replacing by UNK");
              data.get(t)[i] = Dataset.unknownCode;
            }
          }
          else
            data.get(t)[i] = Dataset.unknownCode;
        i++;
      }
//        else
//          throw new RuntimeException("For now, only gap-free alignments accepted");
      final int nSites = align.columns().size();
      final int repeatN=repeats;
      final int nChars = indexer.size()*repeatN;

      final ObservationDimensions od = new ObservationDimensions() {
        @Override
        public int nSites()
        {
         return nSites;
        }
        @Override
        public int nCharacter(int site)
        {
          return nChars;
        }
      };
    int nChars0=indexer.size();
    Map<Taxon,double[][]> converted0 = convert(data, od);
  	final Map<Taxon,double[][]> converted=new HashMap<Taxon,double[][]>();
  	for(Taxon tax:converted0.keySet()){
  	   double[][] mat=converted0.get(tax);  	     	   
  	   double[][] matRepeatN=new double[od.nSites()][nChars];
  		   for(int s=0;s<od.nSites();s++)
  			   for(int c=0;c<nChars0;c++)
  				   for(int l=0;l<repeatN;l++)
     		        matRepeatN[s][c+l*nChars0]=mat[s][c];        		   
  		 converted.put(tax, matRepeatN); 
  	}        	
      
      
      return new Dataset() {
        @Override
        public Map<Taxon, String> getReferenceClusters() { throw new RuntimeException(); }
        @Override
        public boolean hasReferenceClusters()
        {
          return false;
        }
        @Override
        public Map<Taxon, double[][]> observations()
        {
          return converted;
        }
        @Override
        public int nCharacter(int site)
        {
          return nChars;
        }
        @Override
        public int nSites()
        {
          return nSites;
        }
      };
    }
    
    
    private static double[][] aMatrix(Dataset ds) { return ds.observations().values().iterator().next(); }
    public static int nSites(Dataset ds) { return aMatrix(ds).length; }
    public static int nCharacters(Dataset ds, int site) { return aMatrix(ds)[site].length; }
    public static Map<Taxon,double[][]> convert(Map<Taxon,int[]> obs, ObservationDimensions dims)
    {
      Map<Taxon,double[][]> result = new HashMap<Taxon,double[][]>();
      for (Taxon lang : obs.keySet())
        result.put(lang, convert(obs.get(lang), dims));
      return result;
    }
    private static double[][] convert(int[] model,
        ObservationDimensions dims)
    {
      double[][] processed = createObsArray(dims);
      for (int s = 0; s < dims.nSites(); s++)
        for (int c = 0; c < dims.nCharacter(s); c++)
          if (model[s] == unknownCode || model[s] == c)
            processed[s][c] = 1.0;
      return processed;
    }
    public static double[][] createObsArray(ObservationDimensions dims)
    {
      double [][] result = new double[dims.nSites()][];
      for (int i = 0; i < result.length; i++)
        result[i] = new double[dims.nCharacter(i)];
      return result;
    }
    public static double[][] log(double[][] ori)
    {
      double [][] processed = new double[ori.length][];
      for (int s = 0; s < processed.length; s++)
      {
        processed[s] = new double[ori[s].length];
        for (int c = 0; c < processed[s].length; c++)
          processed[s][c] = Math.log(ori[s][c]);
      }
      return processed;
    }
//    public Map<Language,double[]> getSiteObservation(int site)
//    {
//      if (!isFinalState()) throw new RuntimeException();
//      Map<Language,double[]> result = new HashMap<Language,double[]>();
//      for (Arbre<CoalescentNode> node : roots.get(0).nodes())
//        if (node.isLeaf())
//        {
//          DiscreteModelCalculator dmc = (DiscreteModelCalculator) (node.getContents().likelihoodModelCache);
//          if (!dmc.isMissing(site))
//          {
//            double [] cacheCopy = dmc.getCacheCopy(site);
//            NumUtils.expNormalize(cacheCopy);
//            // NB: nodeId not null b/c it's a leaf
//            result.put(node.getContents().nodeIdentifier, cacheCopy);
//          }
//        }
//      return result;
//    }

//    public static Map<Language,double[][]> log(Map<Language,double[][]> ori)
//    {
//      Map<Language,double[][]> result = new HashMap<Language,double[][]>();
//      for (Language lang : ori.keySet())
//      {
//        double [][] model = ori.get(lang);
//        double [][] processed = new double[model.length][];
//        result.put(lang,processed);
//        for (int s = 0; s < processed.length; s++)
//        {
//          processed[s] = new double[model[s].length];
//          for (int c = 0; c < processed[s].length; c++)
//            processed[s][c] = Math.log(model[s][c]);
//        }
//      }
//      return result;
//    }
  }
}
