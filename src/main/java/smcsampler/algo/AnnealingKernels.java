package smcsampler.algo;

import bayonet.distributions.Random;

public interface AnnealingKernels<P extends AnnealedParticle>
{
  /**
   * These are assumed to be exact samples. E.g. prior or a unique pseudo-state.
   */
  P sampleInitial(Random random);
  
  /**
   * Sample from a pi_t invariant kernel where t = temperature
   */
  P sampleNext(Random random, P current, double temperature);
  
  /**
   * Whether sampleNext(.) changes the current state in place.
   */
  boolean inPlace();
  
  /**
   * Only requires is proposal is in place.
   */
  P deepCopy(P particle);
}