/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef ISOLATOR_SEQUENCING_BIAS_HPP
#define ISOLATOR_SEQUENCING_BIAS_HPP

#include <boost/random.hpp>
#include <string>

#include "../pos_table.hpp"
#include "common.hpp"
#include "motif.hpp"
#include "samtools/faidx.h"


/** A representation of tabulated sequence bias. */
struct SeqbiasTabulation
{
    SeqbiasTabulation(unsigned int order=1)
        : order(order), offset(0), bias(NULL)
    {

    }

    ~SeqbiasTabulation()
    {
        delete bias;
    }

    unsigned int order;
    pos_t offset;
    kmer_matrix* bias;
    std::vector<double> divergence;
};


/** A representation of sequencing bias for a particular dataset.
 */
class sequencing_bias
{
    public:
        /* Train a new model.
         *
         * Args:
         *   ref_fn: File name of an indexed FASTA file.
         *   T: Positions of reads
         *   max_reads:  How many reads (at most) to use.
         *   L: How many left positions to consider.
         *   R: How many right positions to consider.
         *   task_name: Logger task_name
         *   tabulation: If non-null, tabulate bias and store in this structure.
         *               The 'order' field of this structure should be
         *               initialized, e.g. to 1.
         *   complexity_penalty: Force sparser models by making this number larger.
         */
        sequencing_bias(const char* ref_fn,
                        PosTable& T,
                        size_t max_reads,
                        pos_t L, pos_t R,
                        const char* task_name,
                        SeqbiasTabulation* tabulation,
                        double complexity_penalty = 1.0);

        ~sequencing_bias();

        double get_bias(const twobitseq& seq, pos_t pos) const;

        /* Return a string of the model graph in dot format. */
        std::string model_graph() const;

        pos_t getL() const { return L; }
        pos_t getR() const { return R; }

    private:
        sequencing_bias();

        void clear();

        void build(const char* ref_fn,
                   PosTable& T,
                   size_t max_reads,
                   pos_t L, pos_t R,
                   const char* task_name,
                   SeqbiasTabulation* tabulation,
                   double complexity_penalty);

        /* left and right sequence context */
        pos_t L, R;

        /* reference sequence */
        faidx_t*    ref_f;
        std::string ref_fn;

        /* trained motif */
        motif* M;

        /* random number generator */
        rng_t rng;
        boost::random::uniform_int_distribution<pos_t> random_uniform_int;
};


#endif
