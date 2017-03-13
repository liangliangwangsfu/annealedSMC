Summary
-------

<!-- [![Build Status](https://travis-ci.org/alexandrebouchard/phylosmcsampler.png?branch=master)](https://travis-ci.org/alexandrebouchard/phylosmcsampler) -->

phylosmcsampler is ...

phylosmcsampler stands for ...


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
  compile group: 'ca.sfu.stat', name: 'phylosmcsampler', version: '1.0.0'
}
```

### Compile using the provided gradle script

- Check out the source ``git clone git@github.com:liangliangwangsfu/phylosmcsampler.git``
- Compile using ``./gradlew installDist``
- Add the jars in ``build/install/phylosmcsampler/lib/`` into your classpath

### Use in eclipse

- Check out the source ``git clone git@github.com:liangliangwangsfu/phylosmcsampler.git``
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root
  - Deselect ``Copy projects into workspace`` to avoid having duplicates


Usage
-----

### Quick start

...
