package pty.smc;

import java.io.*;
import java.util.*;

import pty.smc.PartialCoalescentState;
import nuts.util.CollUtils.*;
import nuts.util.Counter;
import static nuts.util.CollUtils.*;
import static nuts.io.IO.*;
import static nuts.util.MathUtils.*;


public  class LazyPS
{
  private PartialCoalescentState state;
  private double queuedLeftBranch = -1, queuedRightBranch = -1;
  private int queuedLeft = -1, queueRight = -1;
  public LazyPS(PartialCoalescentState pcs) { this.state = pcs; }
  public LazyPS(PartialCoalescentState pcs,  int left, int right, double leftIncrement, double rightIncrement) 
  { 
    this.queuedLeft = left;
    this.queueRight = right;
    this.queuedLeftBranch = leftIncrement;
    this.queuedRightBranch=rightIncrement;
    this.state = pcs; 
  }
  public int nIterationsLeft()
  {
    int result = state.nIterationsLeft();
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
      state = state.coalesce(queuedLeft, queueRight,0.0, queuedLeftBranch, queuedRightBranch);
      queuedLeft = -1;
    }
  }
  public PartialCoalescentState getState()
  {
    applyQueue();
    return state;
  }
}