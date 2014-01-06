

```
     ___           _       _
    |_ _|___  ___ | | __ _| |_ ___  _ __
     | |/ __|/ _ \| |/ _` | __/ _ \| '__|
     | |\__ \ (_) | | (_| | || (_) | |
    |___|___/\___/|_|\__,_|\__\___/|_|
```

Isolator quickly and accurately estimates isoform abundance in RNA-Seq
experiments with MCMC.




## Transcripts Per Truncated Million (TPTM)

RNA-Seq measures relative, not absolute expression. Common measurements of
abundance like FPKM (fragments per exonic-kilobase per million mapped reads) and
TPM (transcripts per million) capture the compositional nature of the data.

A goal of RNA-Seq experiments is often to detect a change in the abundance of a
transcript or gene. This is hampered since a a change in relative abundance does
not necessarily correspond to a change in absolute abundance. As an analogy, if
Bill Gates suddenly got richer, we would all be relatively poorer without any
change in amounts on our bank statements. In any heavy-tailed distribution,
wether it be the distribution of wealth or the distribution of gene expression,
small changes in in tails can dramatically alter relative or compositional
measurements.

Because of this a number of robust normalization procedures were developed for
RNA-Seq accounting for the heavy tails. These procedures seem to work quite
well, but have created a situation where analysis software still reports values
as RPKM, FPKM, or TPM 




