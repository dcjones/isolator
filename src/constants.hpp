
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
}

#endif

