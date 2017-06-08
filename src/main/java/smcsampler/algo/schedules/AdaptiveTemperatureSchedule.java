package smcsampler.algo.schedules;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;

import bayonet.smc.ParticlePopulation;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import smcsampler.algo.Particle;
import smcsampler.algo.SMCStaticUtils;

public class AdaptiveTemperatureSchedule implements TemperatureSchedule
{
  @Arg             @DefaultValue("true")
  public boolean useConditional = true;
  
  @Arg       @DefaultValue("0.999")
  public double threshold = 0.999;
  
  @Override
  public double nextTemperature(ParticlePopulation<? extends Particle> population, double temperature)
  {
    UnivariateFunction objective = objective(population, temperature);
    if (objective.value(1.0) >= 0)
      return 1.0;
    PegasusSolver solver = new PegasusSolver();
    return solver.solve(100, objective, temperature, 1.0);
  }

  private UnivariateFunction objective(ParticlePopulation<? extends Particle> population, double temperature)
  {
    double previousRelativeESS = useConditional ? Double.NaN : population.getRelativeESS();
    return useConditional ? 
        (double proposedNextTemperature) -> SMCStaticUtils.relativeESS(population, temperature, proposedNextTemperature, true)  - threshold:
        (double proposedNextTemperature) -> SMCStaticUtils.relativeESS(population, temperature, proposedNextTemperature, false) - threshold * previousRelativeESS;
  }
}