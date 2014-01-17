

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

  * It implements a full hierarchical Bayesian model of an entire RNA-Seq
    *experiment*. Pooling information leads more accurate expression estimates,
    particularily at low-levels. Conditions can be compared or clustered
    without throwing anything out.
  * It asks the right question, estimating posterior probabilities of effect
    sizes, rather than assigning p-values to a vaguely defined notion of
    differential expression.
  * It models, and is capable of distinguishing between, changes in
    *transcription* and changes in *splicing* between two or more conditions.
  * It models and corrects for sequence bias and GC-bias.
  * Compared to other MCMC approaches, it is exceedingly efficient.

*Note*: I've run this tool on a lot of data, evaluated it with multiple
benchmarks, and am writing a paper about it, but it's still alpha software.
Contact me if you'd like help using it to analyze your data.

# Motivation

Cells use a variety of methods to regulate the expression of RNA products:
modulating transcription, splicing, and degradation, for example. Gene
expression experiments, whether with microarray or sequencing, tend to ignore
the fact that there are multiple possible explanations for a change in the
expression of a transcript. Isolator is an attempt to partially address this
deficiency by building a probabilistic model that is informed by the underlying
biological processes driving expression.

In Isolator, rates of transcription and splicing are included explicitly in a
model that encompasses an entire experiment. Modeling the entire experiment
allows us to pool information and estimate the variability transcription and
splicing within conditions, enabling us to detect when a there is a consistent
change between conditions that is not easy explained by natural variability.

Despite the sophistication of the methods, Isolator is accessible and easy to
use by anyone with a passing familiarity with the command line.

# Installing

Currently, Isolator must be installed by cloning the git repository and building
the source code. It has been tested on OS X and Linux.

## Installing Dependencies

To build Isolator from source, you will need to first install
[Boost](http://www.boost.org/), [HDF5](http://www.hdfgroup.org/HDF5/), and
[cmake](http://www.cmake.org/).

### Linux

On Debian, Ubuntu, or similar.

```sh
apt-get install libboost-dev libhdf5-dev cmake
```

### OSX

The easiest way to install dependencies is with [Homebrew](http://brew.sh/).

```sh
brew install boost hdf5 cmake
```

## Building Isolator

### From git

```sh
git clone https://github.com/dcjones/isolator.git
cd isolator/build
cmake ..
make
make install
```

This will install the `isolator` program to the default destination (usually
`/usr/local/bin`).


# Analyzing an RNA-Seq Experiment

Isolator works in two steps: *analysis* and *summarization*. The analysis step
will run the sampler for a number of iterations collecting estimates of the
model parameters, which it will output to a file. To generate human readable
tables from the raw sampler output,

## What you'll need

Input to isolator is a ![GTF](http://en.wikipedia.org/wiki/Gene_transfer_format)
giving transcript/gene annotations along with one or more
![BAM/SAM](http://en.wikipedia.org/wiki/SAMtools) files giving the aligned reads
for each sample or replicate in the experiment.

Optionally, the genome sequence can also be provided in
![FASTA](http://en.wikipedia.org/wiki/Fasta) format, allowing isolator to
correct for various forms of sequence bias, typically resulting in more accurate
quantification.

## Analyzing the Experiment

Experiments consist of one or more BAM (or SAM) files grouped into conditions.
Isolator puts no restriction on the number of conditions or the number of
samples within conditions, though if the goal is to compare two or more
biological conditions, using multiple replicates is strongy recommended.

The analyze command is invoked like:

```sh
isolator analyze gene_annotations.gtf \
    a1.bam,a2.bam,a3.bam b1.bam,b2.bam,b3.bam
```

BAM files within the same condition are grouped with commas, and conditions are
separated with spaced. Analyze is a general purpose command. It can also be
invoked on a single bam file.

```sh
isolator analyze gene_annotations.gtf a.bam
```

The analyze command run for a while and output a file named (by default)
`isolator-output.h5`. It is an ![HDF5](http://www.hdfgroup.org/HDF5/), which is
a standardized data format than can be accessed using wide variety of tools.
However, data stored in this file is a raw form, and typically far more
information than is needed. To quickly generate to-the-point results, there is
the `isolator summarize` command.


## Summarizing the Analysis

The `isolator analyze` estimates parameters for a full probabilistic model of
the experiment. It is not a simple statistical test that attempts to answer a
single question (e.g. whether a gene is differentialy expressed). As a result, the
results generated from `isolator analyze` can be used to examine the data in
many ways.

This is where the `isolator summarize` comes in. It provides a number of
sub-commands to generate simple tables of results from the sampler output.
What follows are several examples. To see a list of the summarize sub-commands,
run `isolator summarize --list`.

### Transcript and Gene Expression

The simplest sort of summarization one might perform is to generate estimates of
transcript expression.

```sh
isolator summarize transcript-expression isolator-output.h5
```

This will generate a tab-delimited file (named `transcript-expression.tsv`) that looks like:
```csv
gene_name       gene_id transcript_id           sample1_adjusted_tpm    sample2_adjusted_tpm    sample3_adjusted_tpm    sample4_adjusted_tpm    sample5_adjusted_tpm    sample6_adjusted_tpm
mt-Tf   ENSMUSG00000064336      ENSMUST00000082387      3.064260e+01    3.542358e+01    3.068786e+01    3.283976e+01    3.785287e+01    3.647695e+01
mt-Rnr1 ENSMUSG00000064337      ENSMUST00000082388      1.812662e+03    2.566914e+03    2.470344e+03    2.728776e+03    2.445997e+03    2.241069e+03
mt-Tv   ENSMUSG00000064338      ENSMUST00000082389      3.002667e+01    3.085850e+01    3.081490e+01    3.373214e+01    3.507645e+01    3.934946e+01
mt-Rnr2 ENSMUSG00000064339      ENSMUST00000082390      1.337085e+03    1.855240e+03    2.015555e+03    1.996198e+03    1.914501e+03    1.927667e+03
```

Each column gives a point estimate of the adjusted transcripts per million
(TPM). These numbers are adjusted using upper-quantile normalization to make
them comparable across samples.

Credible intervals (i.e. Bayesian confidence interval) for each point estimate
can also be generated by passing the `--credible` argument. The following
command will generate a table that include lower and upper bounds of 95%
credible interval for each point estimate.
```sh
isolator summarize transcript-expression --credible=0.95 isolator-output.h5
```

Expression can also be summarized on the gene level:
```sh
isolator summarize transcript-expression isolator-output.h5
```

### Differential Expression, Transcription, and Splicing

Since Isolator follows a Bayesian methodology, it does not assign p-values to
differential expression, but rather estimates a probability that that
the condition-wise tendendency has changed by more than some minimum fold
change. The researcher may choose the minimum fold change (it is set to a log2 fold
change of 1.5 by default).

The command to generate a table of differential expression results:
```sh
isolator summarize differential-transcript-expression isolator-output.h5
# OR
isolator summarize differential-gene-expression isolator-output.h5
```

The minmum effect size can be set with the `--effect-size` argument.
E.g.  passing `--effect-size 2` will cause the summarize command to compute the
probability of a 2-fold or greater change in expression between conditions.

Expression is a general term for the relative number of copies of an mRNA. A
change in expression can be driven by a change in transcription or in
post-transcriptional processing (e.g. splicing). Transcriptional and
post-transcriptional effects are explicitly modeled in Isolator and can be
separately tested.

To summarize condition pairwise changes transcription:
```sh
isolator summarize differential-transcription isolator-output.h5
```

to summarize condition pairwise changes in splicing:
```sh
isolator summarize differential-splicing isolator-output.h5
```

These commands similarly accept the `-a`, `-b`, and `--effect-size` arguments to
set the conditions being tested and the minimum effect size, respectively.

Note: for simplicitly the term "splicing" is in used here, when in fact it will
capture any change in expression dynamics occouring after transcription
initiation, including changes in splicing and transcription termination.


