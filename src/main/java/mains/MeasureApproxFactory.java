package mains;

import java.io.Writer;
import java.util.Optional;

import baselines.MeasureApproximation;
import bayonet.smc.ParticlePopulation;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import blang.inits.experiments.Experiment;
import briefj.BriefIO;
import smcsampler.algo.AnnealedSMC;
import smcsampler.algo.Particle;


public class MeasureApproxFactory<P extends Particle> extends Experiment
{
  @Arg                                         @DefaultValue("AnnealedSMC")
  public MeasureApproximation<P> approximationAlgorithm = new AnnealedSMC<>();
  
  @SuppressWarnings("rawtypes")
  @Arg         @DefaultValue("HMM")
  public Model model = new HMM();
  
  @Arg                 @DefaultValue("false")
  public boolean forbidOutputFiles = false;
  
  static class EvaluationContext
  {
    final Writer out;
    public EvaluationContext(Writer out)
    {
      this.out = out;
      BriefIO.println(out, "statistic,approximation,truth");
    }
    public void register(String name, double approximation, Optional<Double> truth)
    {
      BriefIO.println(out, name + "," + approximation + "," + (truth.isPresent() ? truth.get() : "NA"));
    }
  }
  
  @Override
  public void run()
  {
    buildAndRun();
  }
  
  @SuppressWarnings("unchecked")
  public ParticlePopulation<P> buildAndRun()
  {
    approximationAlgorithm.setKernels(model.kernels());
    long wallClockTime = System.currentTimeMillis();
    ParticlePopulation<P> approximation = approximationAlgorithm.getApproximation();
    wallClockTime = System.currentTimeMillis() - wallClockTime;
    if (!forbidOutputFiles)
    {
      BriefIO.write(results.getFileInResultFolder("timing.tsv"),
          "wallClockTimeMillis\t" + wallClockTime);
      model.evaluate(approximation, new EvaluationContext(results.getAutoClosedBufferedWriter("results.csv")));
    }
    return approximation;
  }
  
  public static void main(String [] args)
  {
    Experiment.startAutoExit(args); 
  }
}
