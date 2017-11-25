package smcsampler;

import java.util.List;
import java.util.Random;

import fig.basic.Pair;

/**
 * The minimum behavior that users who want to use SMC Sampler should implement
 * 
 * Note: when implementing inference algorithms, one should try to accommodate not only
 * SMCSamplerKernel but also, after cast check, the behavior of the super interface of this interface
 * 
 * @param <S>
 */
public interface SMCSamplerKernel<S>
{
  /**
   * returns null if there are no next state possible (e.g. because of a
   * structured base measure)
   * @param rand
   * @param current
   * @return (state,LOG weight update) of the next state sampled
   * 
   * WARNING: this should be thread-safe: this will be called by several different thread in distributed SMC implementation
   */
  public Pair<S,Double> next(Random rand, S current);
  /**
   * @param partialState
   * @return number of iteration left before partialState becomes a full state
   */
  //public int nIterationsLeft(S partialState);
  /**
   * Should be called only once at the beginnning of the inference algorithm
   * @return
   */
  public S getInitial();
  public  void setCurrentIter(int  currentIter);
  public  int getCurrentIter();
  public void setTemperatureDifference(double temperatureDifference);
  //public void setDeterministicTemperatureDifference(List<Double> temperatureDifference);
  //public double getDeterministicTemperatureDifference(int t);
  public boolean isLastIter();
  public double getDefaultTemperatureDifference(); 
  public double getTemperatureDifference();
  public double getTemperature();
  public void setTemperature(double temperature); 
}