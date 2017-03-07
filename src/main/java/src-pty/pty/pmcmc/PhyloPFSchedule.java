package pty.pmcmc;

import monaco.process.ProcessSchedule;
import monaco.process.ProcessScheduleContext;

public class PhyloPFSchedule implements ProcessSchedule{
	  
	  @Override
	  public void monitor(ProcessScheduleContext context)
	  {
	  }
	  

	  @Override
	  public boolean shouldProcess(ProcessScheduleContext context)
	  {
	    return context.isLastIteration();   
	  }

}
