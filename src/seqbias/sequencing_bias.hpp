/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef ISOLATOR_SEQUENCING_BIAS_HPP
#define ISOLATOR_SEQUENCING_BIAS_HPP

#include <string>

#include "../pos_table.hpp"
#include "common.hpp"
#include "motif.hpp"
#include "samtools/faidx.h"
#include <gsl/gsl_rng.h>

/** A representation of sequencing bias for a particular dataset.
 */
class sequencing_bias
{
    public:
        /* Train a new model.
         *
         * \param ref_fn    File name of an indexed FASTA file.
         * \param T         Positions of reads
         * \param max_reads How many reads (at most) to use.
         * \param L         How many left positions to consider.
         * \param R         How many right positions to consider.
         * \param complexity_penalty Force sparser models by making this number
         *                           larger.
         */
        sequencing_bias(const char* ref_fn,
                        PosTable& T,
                        size_t max_reads,
                        pos_t L, pos_t R,
                        const char* task_name,
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
                   double complexity_penalty);

        /* left and right sequence context */
        pos_t L, R;

        /* reference sequence */
        faidx_t*    ref_f;
        std::string ref_fn;

        /* trained motif */
        motif* M;

        /* random number generator */
        gsl_rng* rng;
};


/** Tabulate the bias in a given dataset, optionally correcting for bias using
 *  the given model. */
kmer_matrix tabulate_bias(double* kl,
                          pos_t L, pos_t R, int k,
                          const char* ref_fn,
                          const char* reads_fn,
                          bool mate1 = true,
                          bool mate2 = true,
                          const char* model_fn = NULL);


#endif
