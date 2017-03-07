package pty.smc;
import java.io.*;
import java.util.*;

import javax.management.RuntimeErrorException;

import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public  class LazyPCS
{
  private PartialCoalescentState state = null, statecopy=null;
  private double queuedMerger = 0;
  private int queuedLeft = -1, queueRight = -1;
  private double queuedDeltaLeft = 0, queuedDeltaRight = 0;
  public LazyPCS(PartialCoalescentState pcs) { this.statecopy = pcs; }
  public LazyPCS(PartialCoalescentState pcs, double delta, int left, int right) 
  { 
    this.queuedLeft = left;
    this.queueRight = right;
    this.queuedMerger = delta;
    this.state = pcs; 
  }
  public LazyPCS(PartialCoalescentState pcs, double delta, int left, int right, double deltaLeft, double deltaRight) 
  {	  
    this.queuedLeft = left;
    this.queueRight = right;
    this.queuedMerger = delta;
    this.state = pcs; 
    this.queuedDeltaLeft = deltaLeft;
    this.queuedDeltaRight = deltaRight;
  }
  public int nIterationsLeft()
  {
    synchronized (this)
    {
    	if (state != null && statecopy != null)
    		throw new RuntimeException();
    	
    	if (statecopy != null)
    		return statecopy.nIterationsLeft();
    	else
    		return state.nIterationsLeft() - 1;
//      int result = state.nIterationsLeft();
//      if (queuedLeft != -1)
//        result--;
////      System.out.println(result);
//      return result;
    }
  }
  public void applyQueue()
  {
    if (queuedLeft == -1) return;	  
    synchronized (this)
    {
    	if (queuedLeft == -1) return;    	
      statecopy = state.coalesce(queuedLeft, queueRight, queuedMerger, queuedDeltaLeft, queuedDeltaRight);
      state=null;
      queuedLeft = -1;
    }
  }
  public PartialCoalescentState getState()
  {
    applyQueue();
    if (statecopy == null)
    	throw new RuntimeException();
    return statecopy;
  }

//public PartialCoalescentState getState()
//{
//    synchronized (this)
//    {
//      if (queuedLeft == -1) return state;    	
//      else{
//    	  
//        state=state.coalesce(queuedLeft, queueRight, queuedMerger, queuedDeltaLeft, queuedDeltaRight);
//        queuedLeft = -1;
//      }
//      return state;
//    }
//}

  
//  public PartialCoalescentState peekState()
//  {
//    synchronized (this)
//    {
//      if (queuedLeft == -1) return state;    	
//      else
//        return state.coalesce(queuedLeft, queueRight, queuedMerger, queuedDeltaLeft, queuedDeltaRight);
//    }
//  }
}