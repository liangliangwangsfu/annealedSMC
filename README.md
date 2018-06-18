Summary
-------

<!-- [![Build Status](https://travis-ci.org/alexandrebouchard/phylosmcsampler.png?branch=master)](https://travis-ci.org/alexandrebouchard/phylosmcsampler) -->

AnnealedSMC is ...

AnnealedSMC stands for ...


Installation
------------


There are several options available to install the package:

### Integrate to a gradle script

Simply add the following lines (replacing 1.0.0 by the current version (see git tags)):

```groovy
repositories {
 mavenCentral()
 jcenter()
 maven {
    url "http://people.stat.sfu.ca/~lwa68/maven/"
  }
}

dependencies {
  compile group: 'ca.sfu.stat', name: 'annealedSMC', version: '1.0.0'
}
```

### Compile using the provided gradle script

- Check out the source ``git clone git@github.com:liangliangwangsfu/annealedSMC.git``
- Compile using ``./gradlew installDist``
- Add the jars in ``build/install/annealedSMC/lib/`` into your classpath

### Use in eclipse

- Check out the source ``git clone git@github.com:liangliangwangsfu/annealedSMC.git``
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root
  - Deselect ``Copy projects into workspace`` to avoid having duplicates


Usage
-----

### Quick start

-useDataGenerator true   -nThousandIters 0.001  -nTax  15 -len  200 -sequenceType DNA -generateDNAdata true -useDataGen4GTRGammaI false -nThreads 2  -treeRate 10 -deltaProposalRate 10 -useNonclock true -useSlightNonclock false -sdScale 0.3 
-iterScalings  100   -methods   ANNEALING   -resamplingStrategy ESS  -nAnnealing 10000 -nSubsampling 10000 -alphaSMCSampler 0.9999    -nSitesPerIndex  10
-essRatioThreshold 0.5 -adaptiveTempDiff true  -runDSMCusingadaptiveTemp  false
-adaptiveType 0     -csmc_trans2tranv 2.0   -mb_trans2tranv 2.0 -setJC false  -fixtratioInMb true  -treePrior unconstrained:exp(10)     -fixNucleotideFreq true   -nReplica  1   -repPerDataPt   10  -mainRand 452  -gen.rand 345 -useCESS true -useNNI true  -useLIS  false  -usenewSS false  -useRevSS false -ntempSS  50  -mcmcfac  1 -usenewDSMC  true
-mrBayesPath  /Users/oudomame/Dropbox/phyloSoftware/mrbayes-3.2.6/src//mb  -neighborPath /Users/oudomame/Dropbox/phyloSoftware/phylip-3.69/exe//neighbor

For settings in our experimental results, please refers to ``https://github.com/shijiaw/annealingExperiment''.
