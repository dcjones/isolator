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

/** A representation of sequencing bias for a particular dataset.
 */
class sequencing_bias
{
    public:
        /** Load a model that has been previously trained, without a reference
         * sequence. (It can not be actually used this way.)
         */
        //sequencing_bias(const char* model_fn);


        /** Load a model that has been previously trained. */
        //sequencing_bias(const char* ref_fn,
                        //const char* model_fn);

        /** Train a new model.
         *
         * \param ref_fn    File name of an indexed FASTA file.
         * \param M1        Positions of mate1 reads.
         * \param M1        Positions of mate2 reads.
         * \param max_reads How many reads (at most) to use.
         * \param L         How many left positions to consider.
         * \param R         How many right positions to consider.
         * \param complexity_penalty Force sparser models by making this number
         *                           larger.
         */
        sequencing_bias(const char* ref_fn,
                        PosTable& T1, PosTable& T2,
                        size_t max_reads,
                        pos_t L, pos_t R,
                        double complexity_penalty = 2.0);


        /** destructor */
        ~sequencing_bias();

        /** Compute the bias across the given region. 
         *
         * The vector returned must be freed with 'delete []'.
         */
        double* get_bias(const char* seqname,
                         pos_t start, pos_t end, strand_t strand) const;

        double get_bias(const twobitseq& seq, pos_t pos) const;

        double* get_mate1_bias(const char* seqname,
                               pos_t start, pos_t end, strand_t strand) const;
        double get_mate1_bias(const twobitseq& seq, pos_t pos) const;

        double* get_mate2_bias(const char* seqname,
                               pos_t start, pos_t end, strand_t strand) const;
        double get_mate2_bias(const twobitseq& seq, pos_t pos) const;


        /** Serialize the model to a file. */
        //void save_to_file(const char* fn) const;

        /** Serialize the model and emit in YAML format. */
        //void to_yaml(YAML::Emitter&) const;

        /** Return a string of the model graph in dot format. */
        std::string model_graph() const;

        pos_t getL() const { return L; }
        pos_t getR() const { return R; }

    private:
        sequencing_bias();

        void clear();

        void build(const char* ref_fn,
                   PosTable& T1, PosTable& T2,
                   size_t max_reads,
                   pos_t L, pos_t R,
                   double complexity_penalty);

        void buildn(motif** Mn,
                    const char* ref_fn,
                    PosTable& T,
                    size_t max_reads,
                    pos_t L, pos_t R,
                    double complexity_penalty,
                    const char* task_name);


        double* get_maten_bias(
                const motif* Mn,
                const char* seqname,
                pos_t start, pos_t end,
                strand_t strand) const;

        double get_maten_bias(
                const motif* Mn,
                const twobitseq& seq, pos_t pos) const;


        /* left and right sequence context */
        pos_t L, R;

        /* reference sequence */
        faidx_t*    ref_f;
        std::string ref_fn;

        /* trained motifs */
        motif* M1; // single-end motif, or double-end mate1 motif
        motif* M2; // double-end mate2 motif

        /* foreground/background distribution over gc content */
        double* gc1[2];
        double* gc2[2];
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
