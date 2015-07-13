
#ifndef ISOLATOR_CONSTANTS_HPP
#define ISOLATOR_CONSTANTS_HPP

#include <cstdlib>

#include "common.hpp"


/* This file contains all of Isolator's dirtly little secrets. All the magic
 * numbers that make it tick but it doesn't like to talk about are layed out
 * here in lurid detail.
 */

namespace constants
{
    /* Number of threads to use. */
    extern unsigned int num_threads;

    /* Number of optimization rounds to produce the maximum posterior estimate.
     * */
    extern unsigned int num_opt_rounds;

    /* During the parameter estimation phase, a work queue of genomic intervals
     * is used to accumulate statistic in parallel. This parameter controls the
     * maximum number of items in that queue. If too large, memory becomes an
     * issue, if too small threads may become starved.
     */
    extern size_t max_estimate_queue_size;

    /* Read this many reads from the BAM file before training seqbias. It will
     * then choose seqbias_num_reads at random from these to train on. */
    extern size_t seqbias_num_collected_reads;

    /* The number of reads to train seqbias with. */
    extern size_t seqbias_num_reads;

    /* The number of positions to the left of the read start to model with
     * seqbias. */
    extern size_t seqbias_left_pos;

    /* The number of positions to the right of the read start to model with
     * seqbias. */
    extern size_t seqbias_right_pos;

    /* Don't use sequences shorter than this to train gc bias */
    extern pos_t gcbias_min_seq_len;

    /* Number of bins to use to model GC bias. */
    extern size_t gcbias_num_bins;

    /* GC upper bounds defining bins for GC content bias. */
    extern float gcbias_bins[];

    /* Don't let the gcbias correction get out of hand to the point of
     * risking over/underflow. */
    extern double gcbias_max_bias;

    /* Minimum and maximum transcript length used to train 3' bias. */
    extern pos_t tpbias_min_tlen, tpbias_max_tlen;

    /* Maximum number of transcripts to use to train 3' bias. */
    extern size_t tpbias_max_transcripts;

    /* Library type; how the mates in paired-end data sets were generated. */
    enum libtype_t {
        LIBTYPE_FR,  //  ======>--------<======
        LIBTYPE_RF,  //  <======--------======>
        LIBTYPE_FF   //  ======>--------======>
    };

    extern libtype_t libtype;

    /* When estimating fragment length distributions, throw out anything larger
     * than this, assuming it's an artifact. */
    extern pos_t max_frag_len;

    /* When estimating fragment length distribution, exclude outlier fragment
     * lengths that have probability less that this number. */
    extern float min_frag_len_pr;

    /* When estimating transcript weights, don't both measuring sequence bias,
     * strand bias etc, for fragments whose length has lower probability than
     * this. */
    extern float transcript_len_min_frag_pr;

    /* Minimum number of paired-end reads needed to estimate the fragment length
     * distribution. */
    extern size_t frag_len_min_pe_reads;

    /* A normal distribution with the following parameters is used as the
     * fragment length distribution when no emperical estimate is available. */
    extern double frag_len_mu;
    extern double frag_len_sd;

    /* Discard fragments with weight lower than this, that might otherwise
     * introduce zeros into the posterior probability. */
    extern float min_frag_weight;

    /* Various epsilons used during numerical optimization of posterior
     * probability. */

    /* During optimization we forbid any transcript from having expression zero
     * as we migth otherwise run into issues with inf/nan values. Instead we use
     * the following number as a lower bound, then round anything near it to 0.0
     * after optimization. */
    extern float zero_eps;

    /* When pesudo transcripts are formed to analyze alternative exons, add
     * flanking exons of this length. */
    extern pos_t alt_exon_flank_length;

    /* To account for imprecision in the transcription start and terminotation
     * sites, we extent the 5p and 3p ends of each transcript. */
    extern pos_t transcript_5p_extension;
    extern pos_t transcript_3p_extension;

    /* It's possible for a transcript to get assigned exceedingly low weight,
     * given the fragment length distribution. Then when a single read lands
     * there, it is assumed that the trascript is very highly expressed. More
     * likely is that the annotation is not correct. (E.g. frequently miRNA are
     * annotated, but pre-miRNA is actually what's being sequenced.) This is a
     * minimum transcript weight that is applied to avoid these situations. */
    extern float min_transcript_weight;

    /* Worker threads for the sampler process blocks of this many components and
     * multireads, respectively, at a time. */
    extern unsigned int sampler_component_block_size;

    /* Normalize to this quantile to compare expression fold change. */
    extern double sample_scaling_quantile;

    /* Before choosing the normalization point, discard all but the n most
     * abundant transcripts. This is to account for very large gene annotation
     * sets containing many transcripts with very low expression, which could
     * bias things. */
    extern size_t sample_scaling_truncation;

    /* Avoid NaNs by rounding zero fragment probabilities up to this tiny
     * positive number.  */
    extern float frag_prob_epsilon;

    /* Round miniscule expression values up to this value, since at that
     * scale changes are meaningless. */
    extern float min_expr;

    /* Hyperparameters for the gamma prior on experiment_mean values. */
    extern double analyze_experiment_mean0;
    extern double analyze_experiment_shape0;
    extern double analyze_experiment_shape;

    /* Hyperparameters for the normal prior on experiment_splice_mu values. */
    extern double analyze_experiment_splice_mu0;
    extern double analyze_experiment_splice_sigma0;

    /* Student-t nu parameter for the experiment tgroup and splice
     * distributions. */
    extern double analyze_experiment_splice_nu;

    /* When there is no data, group level variance parameters can converge
     * to near zero, causing the sampler to get stuck. To avoid this, we force
     * variance parameters to be above some small value. */
    extern double analyze_min_splice_sigma;

    /* How much an alignment is allowed to overhang an intron and still be
     * considered valid. This handles the fairly common case of a read that
     * should be spliced instead aligning into the intron by one or two bases
     * which should span the intron.  */
    extern pos_t max_intron_overlap;

    /* Probabilities controlling how the probability of misalignment is
     * computed.
     */
    extern double aligned_mismatch_pr;
    extern double misaligned_mismatch_pr;
    extern double misalign_prior;

    /* Fragbias is computed explicitly for ends of transcripts and
     * interpolated for the middle. */
    extern pos_t fragbias_endlen;

    /* Scale fragbias to avoid tiny numbers. */
    extern float fragbias_scale;

    /* Minimum probability of correct alignment. Alignments below this
     * will be supressed. */
    extern double min_align_pr;

    /* Maximum number of steps using newton's method before we resort to
     * bisection. */
    extern size_t max_newton_steps;
}

#endif

