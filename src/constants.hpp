
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

    /* Minimum allowable map quality. Alignments lower that this are thrown out.
     * */
    extern unsigned int min_map_qual;

    /* During the parameter estimation phase, a work queue of genomic intervals
     * is used to accumulate statistic in parallel. This parameter controls the
     * maximum number of items in that queue. If too large, memory becomes an
     * issue, if too small threads may become starved.
     */
    extern size_t max_estimate_queue_size;

    /* Minimum length of an exonic region to use while estimating the fragment
     * length distribution. */
    extern pos_t min_estimate_exon_length;

    /* Don't collect training reads from exons shorter than this. */
    extern pos_t seqbias_min_exon_length;

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

    /* Seperate seqbias models are trained for the 3' end 5' ends of
     * transcripts. These numbers define the extent of the ends. */
    extern pos_t seqbias_tp_end;
    extern pos_t seqbias_fp_end;

    /* Don't use sequences shorter than this to train gc bias */
    extern pos_t gcbias_min_seq_len;

    /* Number of bins to use to model GC bias. */
    extern size_t gcbias_num_bins;

    /* Don't let the gcbias correction get out of hand to the point of
     * risking over/underflow. */
    extern double gcbias_max_bias;

    /* Minimum and maximum transcript length used to train 3' bias. */
    extern pos_t tpbias_min_tlen, tpbias_max_tlen;

    /* Maximum number of transcripts to use to train 3' bias. */
    extern size_t tpbias_max_transcripts;

    /* Loess smoothing for GC correction. */
    extern double gc_loess_smoothing;

    /* Library type; how the mates in paired-end data sets were generated. */
    enum libtype_t {
        LIBTYPE_FR,  //  ======>--------<======
        LIBTYPE_RF,  //  <======--------======>
        LIBTYPE_FF   //  ======>--------======>
    };

    extern libtype_t libtype;

    /* Blacklist reads with more than this many alignments. */
    extern unsigned int max_alignments;

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

    /* The degree of smoothing to use for the emperical distribution over
     * fragment lengths. */
    extern float frag_len_dist_smoothing;

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

    /* As an extra precation, don't quantify transcripts where more than this
     * percentage of the reads would be rejected by size selection. */
    extern float min_transcript_fraglen_acceptance;

    /* Symmetric-dirchlet prior precision for transcript mixtures. */
    extern float tmix_prior_prec;

    /* Worker threads for the sampler process blocks of this many components and
     * multireads, respectively, at a time. */
    extern unsigned int sampler_component_block_size;
    extern unsigned int sampler_multiread_block_size;

    /* When doing maximum posterior estimation, we halt when either the relative
     * or absolute error are below respective cutoffs. */
    extern float maxpost_rel_error;
    extern float maxpost_abs_error;

    /* Once consecutive iterations of maximum posterior estimation change the
     * overall probability by less than this amount, the process halts. */
    extern float maxpost_abs_peps;

    /* When using newton's method to find slice extents, only allow so many
     * iterations to avoid cycles. */
    extern unsigned int max_newton_iter;

    /* Number of initial burn-in samples to generate and throw away. */
    extern unsigned int sampler_burnin_samples;

    /* Number of hillclimbing iterations. */
    extern unsigned int sampler_hillclimb_samples;

    /* When looking for changes in splicing, a pooled precision prior is
     * conditioned on the number of isoforms in the group. We don't want this
     * conditioning to make things too sparse for groups with many isoforms, so
     * at least this many group with k isoforms must exist to condition on k,
     * otherwise we just group it with k - 1. Yeah, so this is a little hard to
     * explain tersely....
     */
    extern unsigned int min_tss_group_isoforms_conditioning;

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

    /* Hyperparameters for the normal prior on experiment_tgroup_mu values. */
    extern double analyze_experiment_tgroup_mu0;
    extern double analyze_experiment_tgroup_sigma0;

    /* Hyperparameters for the normal prior on experiment_splice_mu values. */
    extern double analyze_experiment_splice_mu0;
    extern double analyze_experiment_splice_sigma0;

    /* Student-t nu parameter for the experiment tgroup and splice
     * distributions. */
    extern double analyze_experiment_tgroup_nu;
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
}

#endif

