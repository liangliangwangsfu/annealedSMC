/**
 * 
 */
package pty.smc.models;
import java.io.Serializable;

import pty.smc.PartialCoalescentState;
import pty.smc.PartialCoalescentState.CoalescentNode;

public interface LikelihoodModelCalculator 
{
  /**
   * This is the method that gets called repeatedly by 
   * PartialCoalescentState when used in the ParticleFilter
   * 
   * should be fast
   * @return
   */
  public double logLikelihood();
  /**
   * When combining two trees in a coalescent event,
   * return a fresh LikelihoodModelCalculator that
   * sits at the given height (delta heights can 
   * be obtained by calling node{1,2}.height
   * @param node1
   * @param node2
   * @param height
   * @param isRoot
   * @return
   */
  public LikelihoodModelCalculator combine(
      LikelihoodModelCalculator node1, 
      LikelihoodModelCalculator node2, 
      double delta1, double delta2, boolean doNotBuildCache);
  
  public double peekCoalescedLogLikelihood(
      LikelihoodModelCalculator node1, 
      LikelihoodModelCalculator node2, 
      double delta1, double delta2);
  /**
   * In some cases, optimizations are possible for reversible process
   * only.  Return false if unsure.
   * @return
   */
  public boolean isReversible();
  public double extendLogLikelihood(double delta);
  
  

}