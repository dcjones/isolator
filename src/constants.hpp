
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

    /* During the parameter estimation phase, a work queue of genomic intervals
     * is used to accumulate statistic in parallel. This parameter controls the
     * maximum number of items in that queue. If too large, memory becomes an
     * issue, if too small threads may become starved.
     */
    extern size_t max_estimate_queue_size;

    /* Minimum length of an exonic region to use while estimating the fragment
     * length distribution. */
    extern pos_t min_estimate_exon_length;

    /* The number of reads to train seqbias with. */
    extern size_t seqbias_num_reads;

    /* The number of positions to the left of the read start to model with
     * seqbias. */
    extern size_t seqbias_left_pos;

    /* The number of positions to the right of the read start to model with
     * seqbias. */
    extern size_t seqbias_right_pos;

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

    /* Mixtures must sum to 1. This is a tolerance giving the absolute allowable
     * divergence from that constraint. */
    extern float simplex_constraint_eps;

    /* Absolute convergence epsilon for the log posterior objective function. */
    extern float max_post_objf_tolerance;

    /* Absolute convergence epsilon for the mixtures */
    extern float max_post_x_tolerance;

    /* Maximum number of objective function evaluations. */
    extern unsigned int max_post_max_eval;

    /* It's possible for a transcript to get assigned exceedingly low weight,
     * given the fragment length distribution. Then when a single read lands
     * there, it is assumed that the trascript is very highly expressed. More
     * likely is that the annotation is not correct. (E.g. frequently miRNA are
     * annotated, but pre-miRNA is actually what's being sequenced.) This is a
     * minimum transcript weight that is applied to avoid these situations. */
    extern float min_transcript_weight;

    /* Symmetric-dirchlet prior precision for transcript mixtures. */
    extern float tmix_prior_prec;
}

#endif

