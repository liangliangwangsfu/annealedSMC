/**
 * 
 */
package pty.smc.models;
import fig.basic.LogInfo;
import fig.basic.NumUtils;
import nuts.maxent.SloppyMath;
import nuts.util.MathUtils;

/**
 * Assumes reversible dmc for now...
 */
public class ForestModelCalculator implements LikelihoodModelCalculator 
{
  private final DiscreteModelCalculator dmc;
  private final double rootHeight, currentHeight;
  private final double languageInventionRate;
  private final double noLangLogLikelihood, withLangLogLikelihood;
  
  public double posteriorNoLanguagePr()
  {
    final double priorWithLanguage = priorWithLanguage(currentHeight);
    final double [] prs = new double[2];
    prs[0] =   noLangLogLikelihood + Math.log(1.0-priorWithLanguage); // no   lang
    prs[1] = withLangLogLikelihood + Math.log(    priorWithLanguage); // with lang
    NumUtils.expNormalize(prs);
    return prs[0];
  }
  
  public static ForestModelCalculator observation(CTMC ctmc, double [][] initCache,
      double rootHeight, double languageInventionRate)
  {    
    LogInfo.logs("Currently assumes reversibility!");
    DiscreteModelCalculator dmc = DiscreteModelCalculator.observation(ctmc, initCache);
    return  new ForestModelCalculator(dmc, rootHeight, 0.0, languageInventionRate, Double.NEGATIVE_INFINITY, 0.0);
  }
  
  public ForestModelCalculator(DiscreteModelCalculator dmc, double rootHeight,
      double currentHeight, double languageInventionRate,
      double noLangLogLikelihood, double withLangLogLikelihood)
  {
    this.dmc = dmc;
    this.rootHeight = rootHeight;
    this.currentHeight = currentHeight;
    this.languageInventionRate = languageInventionRate;
    this.noLangLogLikelihood = noLangLogLikelihood;
    this.withLangLogLikelihood = withLangLogLikelihood;
  }
  
  private double priorWithLanguage(double height)
  {
    if (height >= rootHeight) 
      throw new RuntimeException("Reached height of " + height + " but max is " + rootHeight);
    return 1.0 - Math.exp(-languageInventionRate * (rootHeight - height));
  }

  public double extendLogLikelihood(double delta)
  {
    //final double modifiedWithLangLogLikelihood = dmc.extendLogLikelihood(delta);
    final double prLangInvEvent = 1.0 - Math.exp(-languageInventionRate * delta);
    final double priorWithLanguage = priorWithLanguage(currentHeight + delta);
    return SloppyMath.logAdd(
       // with lang
       Math.log(priorWithLanguage) + withLangLogLikelihood, 
       // without
       Math.log(1.0 - priorWithLanguage) + 
         SloppyMath.logAdd(
             Math.log(prLangInvEvent) + withLangLogLikelihood, // Note: when not reversible, this is actually an approx...
             Math.log(1.0-prLangInvEvent) + noLangLogLikelihood
         )
    );
  }
  
  public Object calculate(
      LikelihoodModelCalculator node1, LikelihoodModelCalculator node2, 
      double v1, double v2,
      boolean isPeek, boolean doNotBuildCache)
  {
    final ForestModelCalculator 
      n1 = (ForestModelCalculator) node1,
      n2 = (ForestModelCalculator) node2;
    MathUtils.checkClose(
        v1 + n1.currentHeight, 
        v2 + n2.currentHeight);
    final double newHeight = v1 + n1.currentHeight;
    double resultNoLang = Double.NEGATIVE_INFINITY;
    final double 
      prLangInvEvent1 = 1.0 - Math.exp(-languageInventionRate * v1),
      prLangInvEvent2 = 1.0 - Math.exp(-languageInventionRate * v2);
    // nolang-nolang
    resultNoLang = SloppyMath.logAdd(resultNoLang, 
        Math.log(1.0-prLangInvEvent1) + n1.noLangLogLikelihood +
        Math.log(1.0-prLangInvEvent2) + n2.noLangLogLikelihood);
    // lang-nolang
    resultNoLang = SloppyMath.logAdd(resultNoLang, 
        Math.log(prLangInvEvent1) + n1.withLangLogLikelihood +
        Math.log(1.0-prLangInvEvent2) + n2.noLangLogLikelihood);    
    // nolang-lang
    resultNoLang = SloppyMath.logAdd(resultNoLang, 
        Math.log(1.0-prLangInvEvent1) + n1.noLangLogLikelihood +
        Math.log(prLangInvEvent2) + n2.withLangLogLikelihood);    
    // lang-lang
    resultNoLang = SloppyMath.logAdd(resultNoLang, 
        Math.log(prLangInvEvent1) + n1.withLangLogLikelihood +
        Math.log(prLangInvEvent2) + n2.withLangLogLikelihood); 
    if (isPeek) 
    {
      final double resultWithLang = dmc.peekCoalescedLogLikelihood(n1.dmc, n2.dmc, v1, v2);
      return logLikelihood(newHeight, resultWithLang, resultNoLang);
    }
    else
    {
      final DiscreteModelCalculator resultDMC = (DiscreteModelCalculator) dmc.combine(n1.dmc, n2.dmc, v1, v2, doNotBuildCache);
      final double resultWithLang = resultDMC.logLikelihood();
      return new ForestModelCalculator(resultDMC, rootHeight, newHeight, languageInventionRate, resultNoLang, resultWithLang);
    }
  }

  private transient double ll = Double.NaN;
  public double logLikelihood()
  {
    if (!Double.isNaN(ll)) return ll;
    ll = logLikelihood(currentHeight, withLangLogLikelihood, noLangLogLikelihood);
    return ll;
  }
  private double logLikelihood(double currentHeight, double withLangLogLikelihood, double noLangLogLikelihood)
  {
    final double priorWithLanguage = priorWithLanguage(currentHeight);
    return SloppyMath.logAdd(
       // with lang
       Math.log(priorWithLanguage) + withLangLogLikelihood,
       // without
       Math.log(1.0 - priorWithLanguage) + noLangLogLikelihood);
  }

  public boolean isReversible()
  {
    return true;
  }

  public double peekCoalescedLogLikelihood(LikelihoodModelCalculator node1,
      LikelihoodModelCalculator node2, double delta1, double delta2)
  {
    return (Double) calculate(node1, node2, delta1, delta2, true, true);
  }

  public LikelihoodModelCalculator combine(
      LikelihoodModelCalculator node1, LikelihoodModelCalculator node2, 
      double v1, double v2, boolean doNotBuildCache)
  {
    return (LikelihoodModelCalculator) calculate(node1, node2, v1, v2, false, doNotBuildCache);
  }

  
}