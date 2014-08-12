/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#ifndef ISOLATOR_MOTIF_HPP
#define ISOLATOR_MOTIF_HPP

#include "twobitseq.hpp"
#include "kmer_matrix.hpp"
#include <deque>

/** A nucleotide motif, modeled with a directed graph (Bayesian network).
 *
 * These can only be trained in pairs: background and foreground.
 */

class motif
{
    public:
        /** Train a new motif, given a set of foreground and background
         * sequences. */
        motif(const std::deque<twobitseq*>& foreground_seqs,
              const std::deque<twobitseq*>& background_seqs,
              size_t m,
              size_t max_parents = 0,
              size_t max_distance = 0,
              double complexity_penalty = 0.5,
              const char* task_name = NULL);

        /** read from yaml */
        //motif(const YAML::Node& node);

        /** copy constructor */
        motif(const motif&);

        ~motif();

        /** Return the likelihood of a sequnce, given the motif. */
        double eval(const twobitseq&, size_t offset) const;

        /** emit yaml serialization */
        //void to_yaml(YAML::Emitter& out) const;

        /** print the model graph in GraphViz (dot) format,
         *  adjusting position labels by the given offset. */
        std::string model_graph(int offset) const;

        /* Estimate the JS divergence between two motifs.
         *
         * They motifs must of the same length. The estimation is performed by
         * sampling over sequences of length 40.
         */
        void estimate_js_divergence(const motif& other, double* d0, double* d1,
                                      size_t num_samples=10000);

    private:
        motif();

        void eval_foreground_background(const twobitseq& seq, size_t offset,
                                        double* p0, double* p1) const;

        size_t num_parents(size_t i) const;
        size_t num_params() const;

        void set_edge(int i, int j, bool x);
        bool has_edge(int i, int j) const;

        size_t m; //< number of positions */

        kmer_matrix* P0; //< background distribution parameters
        kmer_matrix* P1; //< foreground distribution parameters

        /** An n*n 0-1 matrix marking if note i is a parent of j */
        bool* parents;
        size_t* nparents;

        friend class motif_trainer;
};


#endif

