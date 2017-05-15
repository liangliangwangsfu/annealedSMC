package smcsampler;
import static nuts.io.IO.i;
import static nuts.io.IO.writeToDisk;
import static nuts.util.CollUtils.map;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import ev.to.NexusWriter;
import ma.MSAParser;
import ma.MSAPoset;
import ma.RateMatrixLoader;
import ma.SequenceType;
import nuts.io.IO;
import nuts.lang.StringUtils;
import pty.RootedTree;
import pty.RootedTree.RootedTreeProcessor;
import fig.basic.IOUtils;
import fig.basic.Option;


public class MrBayes implements Runnable
{
	public static final MrBayes instance = new MrBayes();

	@Option public String mrBayesPath = "mb";
	@Option public int nMCMCIters = 10000;
	@Option public int nChains = 1;
	@Option public int seed = 1297732343;
	@Option
	public String treePrior = "clock:coalescence"; // unconstrained:exp(10.0)
	@Option public double mbRate = 1.0;
	@Option public boolean setToK2P = false;  
	@Option public boolean setFixCoalescentPr = true;
	@Option public boolean fixNucleotideFreq = false;
	@Option public boolean set2nst = false;
	@Option public boolean setGTRGammaI = false;
	@Option public boolean setInv = false;
	@Option public boolean fixGTRGammaPara=false;
	@Option public double alpha=0.5;
	@Option public double[] stationaryDistribution=new double[]{0.2, 0.26, 0.33, 0.21};
	@Option public double[] subsRates=new double[]{0.26, 0.18, 0.17, 0.15, 0.11, 0.13};
	@Option public double  mb_trans2tranv=2.0; 
	@Option public boolean setstarttree=false;  
	//@Option public boolean setSSinMB=false;
	@Option public boolean fixtratioInMb=false;

	// the next two can also be provided directly as arguments in computeSamples
	@Option public SequenceType st = SequenceType.RNA;
	@Option public File alignmentInputFile = null;

	private File  workingDir = null;
	private String starttreeString="";

	public static final String NEX_FILE = "data.nex";
	public File computeSamples(MSAPoset alignment, SequenceType st)
	{
		workingDir = IO.getTempDir("temp-mrbayes");
		File nexusFile = new File(workingDir, NEX_FILE);
		File mrBayesCmd = new File(workingDir, "mrbayes.cmd");
		NexusWriter.writeNexus(alignment, nexusFile, st);
		writeMrBayesCmd(mrBayesCmd);    
		String msg = IO.call("" + mrBayesPath + " " + mrBayesCmd.getName(), null, workingDir);
//		System.out.println(msg);
		writeToDisk(new File(workingDir, "mrbayes-stdout"), msg);
		return workingDir;
	}

	public String computeMarginalLike(MSAPoset alignment, SequenceType st)
	{
		workingDir = IO.getTempDir("temp-mrbayes-ss");
		File nexusFile = new File(workingDir, NEX_FILE);
		File mrBayesCmd = new File(workingDir, "mrbayes.cmd");
		NexusWriter.writeNexus(alignment, nexusFile, st);
		writeMrBayesCmd_(mrBayesCmd,true);    
		String msg = IO.call("" + mrBayesPath + " " + mrBayesCmd.getName(), null, workingDir);		
		writeToDisk(new File(workingDir, "mrbayes-stdout"), msg);
//		LogInfo.logs("grep Mean: "+mrBayesFolder+"/time=*/mrbayes-stdout |awk {'print $3'}");
		String marginalLikeMean = IO.call("grep Mean: "+workingDir+"/mrbayes-stdout", null, workingDir);		
//		String marginalLikeMean = IO.call("grep 1     -  "+workingDir+"/mrbayes-stdout", null, workingDir);		     
		return marginalLikeMean.substring(marginalLikeMean.indexOf(":")+1, marginalLikeMean.indexOf(".")+3);
	}

	
	public void processMrBayesTrees(RootedTreeProcessor rtp)
	{
		processMrBayesTrees(rtp, 1);
	}

	public void setStartTree(String startTree)
	{	   				
		starttreeString=startTree;
		if(!setstarttree) setstarttree=true;
	}

	public String getStartTree()
	{
		return starttreeString; 
	}

	public void cleanUpMrBayesOutput()
	{
		for (File f : IO.ls(workingDir))
			f.delete();
		workingDir.delete();
	}
	/**
	 * WARNING: MrBayes runs are 1-indexed!
	 * @param rtp
	 * @param runNumber
	 */
	public void processMrBayesTrees(RootedTreeProcessor rtp, int runNumber)
	{		
		processMrBayesTrees(new File(workingDir, "data.nex.run" + runNumber + ".t"), rtp);
	}
	public static void processMrBayesTrees(File file, RootedTreeProcessor rtp)
	{
		Map<String,String> translation = translation(file);
		for (String line : i(file))
			if (line.matches("^\\s*tree rep.*") || line.matches("^\\s*tree gen.*")  )			
				rtp.process(parseTree(line,translation));
	}
	private static Map<String, String> translation(File file)
	{
		Map<String, String> trans = map();
		loop:for (String line : IO.i(file))
			if (newickPattern.matcher(line).matches())
				break loop;
			else if (translPattern.matcher(line).matches())
			{
				List<String> matches = StringUtils.multiSelectFirstRegex(translPattern, line);
				trans.put(matches.get(0), matches.get(1));
			}
		return trans;
	}
	private static final Pattern translPattern = Pattern.compile("\\s*([0-9]+)\\s(.*)[,;]");
    private static final Pattern newickPattern = Pattern.compile(".*[=]\\s(.*)[;]");	
    
	//  private static final Pattern speciesCode = Pattern.compile("[;)]([0-9]+)[:]");
	private static RootedTree parseTree(String line,
			Map<String, String> translation)
	{		
		String newickStr = StringUtils.selectFirstRegex(newickPattern, line);				
		newickStr=newickStr.substring(4); 	
		// hack (cannot translate in tree rep because of limitation of newick parser		
		for (String code : translation.keySet())
		{
			String codeCtx = StringUtils.selectFirstRegex("([,)(])" + code + "[:]", newickStr);			
			newickStr = newickStr.replace(codeCtx + code + ":", codeCtx + translation.get(code) + ":");
		}

		//    for (int i = translation.size(); i >= 1; i--)
		//      newickStr = newickStr.replace("" + i + ":" , translation.get("" + i) + ":"); 

		newickStr = newickStr + ";";		
		return RootedTree.Util.fromNewickString(newickStr);
	}
	
	
	private void writeMrBayesCmd(File mrBayesCmd)
	{
		writeMrBayesCmd_(mrBayesCmd, false);
	}
	
	private void writeMrBayesCmd_(File mrBayesCmd, boolean setSSinMB)
	{
		RateMatrixLoader.DEFAULT_TRANS2TRANV=mb_trans2tranv;    
		final double b4 = RateMatrixLoader.beta(RateMatrixLoader.DEFAULT_TRANS2TRANV) * 4.0;
		final double a4 = RateMatrixLoader.alpha(RateMatrixLoader.DEFAULT_TRANS2TRANV) * 4.0;
		String fixGtrGammaStr=null;
		if(fixGTRGammaPara)
			fixGtrGammaStr="prset statefreqpr=fixed("+stationaryDistribution[0]+","+stationaryDistribution[1]+","+stationaryDistribution[2]+","+stationaryDistribution[3]+");\n"+
					"prset revmatpr=fixed("+subsRates[0]+","+subsRates[1]+","+subsRates[2]+","+subsRates[3]+","+subsRates[4]+","+subsRates[5]+");\n"+
					"prset shape=fixed("+alpha+");\n";
		if(setstarttree)
		{
			File treeCmd = new File(workingDir, "start.tree");
			PrintWriter outtree = IOUtils.openOutEasy(treeCmd);
			outtree.append("execute " + NEX_FILE + ";\n" +"begin trees;\n"+"tree starttree="+starttreeString+"\nend;\n\n");
			outtree.close();
			String cmdStr="cat start.tree | sed 's/internal_[0-9]*_[0-9]*//g'"+  " > mrbayes.cmd";
			IO.call("bash -s",cmdStr,workingDir);      
		}

		PrintWriter out = IOUtils.openOutAppendEasy(mrBayesCmd);
		out.append(
				"begin mrbayes;\n" +
						"set autoclose=yes nowarn=yes;\n" +
						(setstarttree?"":"execute " + NEX_FILE + ";\n") +
						"prset brlenspr=" + treePrior + (mbRate == 1.0  ? "" : "(" + mbRate + ")" ) + ";\n" +
						//        (setFixCoalescentPr ? 
						//          "prset thetapr=fixed(1.0);\n" : "") +
						"mcmcp ngen=" + nMCMCIters + ";\n" +
						"mcmcp Nchains=" + nChains + ";\n" +
//						"mcmcp seed=" + Math.abs(seed) + ";\n" +
						"set seed=" + Math.abs(seed) + ";\n" +  
						"set scientific=no;\n"+
						(setGTRGammaI?" lset nst=6 rates="+(setInv?"invgamma":"gamma")+" ngammacat=4;\n"+"prset tratiopr = beta(1, 1);\n":"")+
						(set2nst?"lset nst=2;\n"+"prset tratiopr = beta(1, 1);\n":"lset nst=6;\n")+
						(fixNucleotideFreq?"prset statefreqpr=fixed(0.25,0.25,0.25,0.25);\n":"")+        
						(setToK2P ? "prset revmatpr=fixed(" + b4 +"," + a4 + "," + b4 + "," + b4 +"," + a4 + "," + b4 + ");\n":"")+
						(fixtratioInMb? "prset tratiopr = fixed("+mb_trans2tranv+");\n":"")+
						(fixGTRGammaPara?fixGtrGammaStr:"")+
						(setstarttree?("startvals  tau = starttree;\n"+
								"startvals  V = starttree;\n"):"")+
						"propset  NNI(Tau,V)$prob=25\n"+
						"propset Multiplier(V)$prob=25\n"+
						"propset Nodeslider(V)$prob=25\n"+
						"propset TLMultiplier(V)$prob=25\n"+
						(setSSinMB?"ss;\n":"mcmc;\n" +								
								"sumt;\n")+								 
				"end;\n");
		out.close();
	}

	public static void main(String [] args)
	{
		//    IO.run(args, new MrBayes());

//		String newickStr="(Selaginella_uncinata:1.000000e-01,((Panax_ginseng:1.000000e-01,Daucus_carota:1.000000e-01):1.000000e-01,Liriodendron_tulipifera:1.000000e-01):1.000000e-01,Coffea_arabica:1.000000e-01);";
		String newickStr="(Selaginella_uncinata:0.1,((Panax_ginseng:0.1,Daucus_carota:0.1):0.1,Liriodendron_tulipifera:0.1):0.1,Coffea_arabica:0.1);";
//		String newickStr="((Liriodendron_tulipifera:0.100000,(Selaginella_uncinata:0.100000,Panax_ginseng:0.100000):0.100000):0.100000,Daucus_carota:0.100000,Coffea_arabica:0.100000);";
		RootedTree rt=RootedTree.Util.fromNewickString(newickStr);
	    System.out.println(rt);
		
		/*
		// debug:
		File mrBayesTrees = new File("/Users/bouchard/w/ptychodus/data/mrBayesOnWalsIE/nAryCharacters-1m-iters/inMrBayesFolder/test.nex.run1.t");
		//    for (String line : i(mrBayesTrees))
		//      System.out.println(line);
		processMrBayesTrees(mrBayesTrees, new RootedTreeProcessor() {

			@Override
			public void process(RootedTree rt)
			{
				System.out.println(rt);
			}
		});
		*/
	}

	@Override
	public void run()
	{
		MSAPoset align = MSAParser.parseMSA(alignmentInputFile);
		computeSamples(align, st);
		processMrBayesTrees(new RootedTreeProcessor() {

			@Override
			public void process(RootedTree rt)
			{
				System.out.println(rt);
			}
		});
	}
}
