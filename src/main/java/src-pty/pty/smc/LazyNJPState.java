package pty.smc;
import java.io.*;
import java.util.*;

import ev.ex.NJPState;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public  class LazyNJPState
{
  private NJPState state;
  private double queuedLeftBranch = -1, queuedRightBranch = -1;
  private int queuedLeft = -1, queueRight = -1;
  public LazyNJPState(NJPState njps) { this.state = njps; }
  public LazyNJPState(NJPState njps, int left, int right, double leftIncrement, double rightIncrement) 
  { 
    this.queuedLeft = left;
    this.queueRight = right;
    this.queuedLeftBranch = leftIncrement;
    this.queuedRightBranch=rightIncrement;
    this.state = njps; 
  }
  public int nIterationsLeft()
  {	  
    int result = state.pcs.nIterationsLeft();
    if (queuedLeft != -1)
      result--;
    return result;
  }
  public void applyQueue()
  {
    if (queuedLeft == -1) return;
    synchronized (this)
    {
      if (queuedLeft == -1) return;
      state = state.coalesce(queuedLeft, queueRight, queuedLeftBranch, queuedRightBranch);
      queuedLeft = -1;
    }
  }
  public NJPState getState()
  {
    applyQueue();
    return state;
  }
}