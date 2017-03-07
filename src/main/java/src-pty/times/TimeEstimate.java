package times;

import java.util.Random;

import fig.basic.LogInfo;
import fig.basic.Option;
import fig.exec.Execution;

public class TimeEstimate implements Runnable {

	@Option public String treeFile = "newick.txt";
	@Option public String timeFile = "time.txt";
	@Option public String rateFile = "rate.txt";
	@Option public String ingroupFile = "ingroup.txt";
	@Option public int maxiters = 10000;
	@Option public int burnin  = 100;
	
	public static final Random rand= new Random(1);
	TreeState ts;

	public static void main (String args []){
	    Execution.monitor = true;
	    Execution.makeThunk = false;
	    Execution.create = true;
	    Execution.useStandardExecPoolDirStrategy = true;
		Execution.run(args, new TimeEstimate()
		);
	}
	
	
		
	
	public void run() {
		// TODO Auto-generated method stub
		ts =  new TreeState (treeFile, timeFile, rateFile, ingroupFile);
		ts.init(rand);
		mcmc mc = new mcmc (ts, rand);
		
		for (int i = 0 ; i < burnin + maxiters; i++ ) {
			TreeState currentState = (TreeState)mc.sample();
			LogInfo.logs("Iteration = " +  i + "\t" + currentState);
		}
	}

}
