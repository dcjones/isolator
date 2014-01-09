

```
     ___           _       _
    |_ _|___  ___ | | __ _| |_ ___  _ __
     | |/ __|/ _ \| |/ _` | __/ _ \| '__|
     | |\__ \ (_) | | (_| | || (_) | |
    |___|___/\___/|_|\__,_|\__\___/|_|
```

Isolator analyzes RNA-Seq experiments.

There are many methods to analyze RNA-Seq data. Isolator differs in few
important ways.

  * It implements a full heirarchical Bayesian model of an RNA-Seq *experiment*.
  * It models and is capable of distinguishing between changes in
    *transcription* and changes in *splicing* between conditions.
  * It models and corrects for sequence bias and GC-bias.
  * Compared to other MCMC approaches, it is exceedingly efficient.


# Motivation

Cells use a variety of methods to regulate the expression of RNA products:
modulating transcription, splicing, and degredation, for example. Gene
expression experiments, whether with microarray or sequencing, tend to ignore
the fact that there are multiple possible explainations for a change in the
expression of a transcript. Isolator is an attempt to partialy address this
deficiency by building a probabalistic model that is informed by the underlying
biological processes driving expression.

In Isolator, rates of transcription and splicing are included explicitly in a
model that encompasses an entire experiment. Modeling the entire experiment
allows us to pool information and estimate the variability transcription and
splicing within conditions, enabling us to detect when a there is a consistent
change between conditions that is not easy explained by natural variability.

Despite the sophistication of the methods, Isolator is accessable and easy to
use by anyone with a passing familiarity with the command line.

# Installing

Currently, Isolator must be installed by cloning the git repository and building
the source code. It has been tested on OS X and Linux.

## Installing Dependencies

To build Isolator from source, you will need to first install
[Boost](http://www.boost.org/), [HDF5](http://www.hdfgroup.org/HDF5/), and
[cmake](http://www.cmake.org/).

### Linux

On Debian or Ubuntu

```sh
```

### OSX

The easiest way to install dependencies is with [Homebrew](http://brew.sh/).

```sh
brew install boost hdf5 cmake
```

## Building Isolator

### From git


### From a tarball


# Analyzing an RNA-Seq Experiment

TODO




