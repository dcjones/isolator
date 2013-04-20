/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#ifndef TWOBITSEQ_MATRIX_HPP
#define TWOBITSEQ_MATRIX_HPP

#include <cstdlib>
#include <string>

#include "common.hpp"
#include "kmer_matrix.hpp"


/** a nucleotide sequence encoded packed two-bits per nucleotide.
 *
 * This is useful for a couple reasons:
 *   1. this uses a fourth of the memory as char arrays would use.
 *   2. we need not worry about capitalization, A vs U, N's, or the like.
 */
class twobitseq
{
    public:
        /** create an empty sequence. */
        twobitseq();

        /** create a sequence from a DNA/RNA sequence string */
        twobitseq(const char* seq);

        /** copy constructor */
        twobitseq(const twobitseq&);

        ~twobitseq();

        /** Return the length of the sequence. */
        size_t size() const;

        /** Shrink or enlarge the sequence.
         *
         * If enlarged, the extra nucleotides will be uninitialized */
        void resize(size_t n);

        /** Compute the reverse complement in place. */
        void revcomp();

        void operator = (const twobitseq&);
        void operator = (const char* seq);

        /* get the kmer ending at position i */
        kmer get_kmer(int k, pos_t i);

        /** extract a kmer.
         *
         * The kmer is made from all the positions with a '1' bit in the 'mask'
         * vector, at the given offset in the sequence.
         */
        int make_kmer(kmer& K, size_t offset, bool* mask, size_t mask_len) const;

        /** Copy a part of another sequence. */
        void copy(const twobitseq& src, pos_t src_start, pos_t dest_start, pos_t n);

        /** Count number of Gs and Cs between i and j inclusively. */
        size_t gc_count(pos_t i, pos_t j) const;

        /* Count number of Gs and Cs in the entire sequence. */
        size_t gc_count() const;

        kmer getnuc(size_t i) const;
        void setnuc(size_t i, kmer K);

        std::string to_string() const;

    private:

        kmer* xs;
        size_t n;

        /** how many nucleotides can we back in a kmer (or, 'k') */
        static const size_t max_kmer;
};


/** Convert a nucleotide charactor to a number, using the same scheme as
 * twobitseq */
kmer nuc_to_num(char c);

/** Convert a number n encoding a kmer into a string of nucleotides */
void num_to_nuc(char* dest, kmer K, int k);

kmer twobitcomp(kmer K);


#endif

