package pty.io;
import goblin.Taxon;
import java.io.File;
import java.util.List;
import java.util.Map;

import nuts.util.CollUtils;
import ma.MSAPoset;

public class ExtractSubDataset { 
	 public static void subset(File alignmentFile, List<Taxon> taxa2remove, File  output)
	    {		 		 
	    MSAPoset align = MSAPoset.parseAlnOrMsfFormats(alignmentFile);
	    Map<Taxon,String> seqMap=align.sequences();
	    Map<Taxon,String> newMap=CollUtils.map(); 
	      for (Taxon t : align.sequences().keySet())
	    	  if(!taxa2remove.contains(t))
	    	  {
	    	   System.out.println(seqMap.get(t));	  
	    		  newMap.put(t, seqMap.get(t));
	    	  }
	    	  
	       MSAPoset msa=new MSAPoset(newMap);    	  
	       
	      msa.toMultiAlignmentObject().saveToMSF(output);        
	    }
	
		public static void main(String [] args)
		{
		List<Taxon> taxa2remove=CollUtils.list();
		//---------------------- clade A
		taxa2remove.add(new Taxon("Metriaclima_zebra"));
		taxa2remove.add(new Taxon("Buccochromis_lepturus"));
		taxa2remove.add(new Taxon("Champsochromis_spilorhynchus"));
		taxa2remove.add(new Taxon("Lethrinops_auritus"));
		taxa2remove.add(new Taxon("Rhamphochromis_sp"));

		//---------------------- clade B
		taxa2remove.add(new Taxon("Lobochilotes_labiatus"));
		taxa2remove.add(new Taxon("Petrochromis_orthognathus"));
		taxa2remove.add(new Taxon("Gnathochromis_pfefferi"));
		taxa2remove.add(new Taxon("Tropheus_moorii"));

		//---------------------- clade C
		taxa2remove.add(new Taxon("Callochromis_macrops"));
		taxa2remove.add(new Taxon("Cardiopharynx_schoutedeni"));
		taxa2remove.add(new Taxon("Ophthalmotilapia_ventralis"));
		taxa2remove.add(new Taxon("Xenotilapia_flavipinnus"));
		taxa2remove.add(new Taxon("Xenotilapia_sima"));

//		//---------------------- clade E
		taxa2remove.add(new Taxon("Perissodus_microlepis_T32a"));
		taxa2remove.add(new Taxon("Perissodus_microlepis_T51a"));
		taxa2remove.add(new Taxon("Cyphotilapia_frontosa"));
		taxa2remove.add(new Taxon("Tanganicodus_irsacae"));
		taxa2remove.add(new Taxon("Limnochromis_auritus"));
		taxa2remove.add(new Taxon("Paracyprichromis_brieni"));
		
//		---------------------- clade F
//		taxa2remove.add(new Taxon("Oreochromis_niloticus"));
//		taxa2remove.add(new Taxon("Tylochromis_polylepis"));
//		taxa2remove.add(new Taxon("Boulengerochromis_microlepis"));
//		taxa2remove.add(new Taxon("Bathybates_sp"));
		taxa2remove.add(new Taxon("Cichlasoma_citrinellum"));
	
		File alignmentFile=new File("/Users/liangliang/data/cichlidFish/cichlidFishWithSpeciesName.msf");
		//File  output=new File("/Users/liangliang/data/cichlidFish/subsetCichlidFishWithSpeciesName.fasta");
		//File alignmentFile=new File("/Users/liangliang/data/cichlidFish/subsetCichlidFish.msf");  //tribes A, C, D. 
		//File  output=new File("/Users/liangliang/data/cichlidFish/subsetCichlidFishTribeAD.msf");
		//File  output=new File("/Users/liangliang/data/cichlidFish/subsetCichlidFishTribeCD.msf");
		//File  output=new File("/Users/liangliang/data/cichlidFish/subsetCichlidFishTribeD.msf");
		//File  output=new File("/Users/liangliang/data/cichlidFish/subsetCichlidFishTribeDEF.msf");
		File  output=new File("/Users/liangliang/data/cichlidFish/subsetCichlidFishTribeDF.msf");
		subset(alignmentFile, taxa2remove, output); 
		}    
		
}
