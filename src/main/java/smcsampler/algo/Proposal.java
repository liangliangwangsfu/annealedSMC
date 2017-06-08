package smcsampler.algo;

import java.util.Random;

public interface Proposal<P extends Particle>
{
  /**
   * These are assumed to be exact samples. E.g. prior or a unique pseudo-state.
   */
  P sampleInitial(Random random);
  P propose(Random random, P current, double temperature);
  
  boolean inPlace();
  /**
   * Only requires is proposal is in place.
   */
  P deepCopy(P particle);
}