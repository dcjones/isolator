/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#ifndef KMER_MATRIX_HPP
#define KMER_MATRIX_HPP

//#include "yaml-cpp/yaml.h"
#include <deque>


typedef unsigned int kmer;

/**
 * A class that stores a matrix over kmers.
 *
 * Columns are indexed by k-mer (encoded in two bits), and there are as many
 * columns as there are k-mers (i.e., 4^k). Any number of rows are allowed, so
 * that the matrix might be visualize as,
 *
 *     AAA AAT AAC AAG ATA ...
 *  1  a11 a12 a13 ...
 *  2  a21 a22 a23 ...
 *  3  a31 a32 a33 ...
 *  ...   ...
 *  n  an1 an2 an3 ...
 * 
 */
class kmer_matrix
{
    public:
        /** read a matrix from a yaml file. */
        //kmer_matrix(const YAML::Node& node);

        /** create an empty matrix with n rows, over all k-mers. */
        kmer_matrix(size_t n, size_t k);

        /** copy constructor. */
        kmer_matrix(const kmer_matrix&);

        ~kmer_matrix();

        /** standard assignment */
        void operator=(const kmer_matrix&);

        /** write tha matrix to a yaml file */
        //void to_yaml(YAML::Emitter& out) const;

        /** get element i, j */
        double& operator()(size_t i, size_t j);

        /** set every element in the matrix to x */
        void set_all(double x);

        /** copy row i into vector x */
        void get_row(size_t i, double* xs);

        /** copy a vector x into row i */
        void set_row(size_t i, double* xs);

        /** set every element in row i to x */
        void set_row(size_t i, double x);

        /** return the number of rows (or, 'n'). */
        size_t nrows() const;

        /** return the number of columns (or, 'm'), which is always a power of 4
         * */
        size_t ncols() const;

        /** return 'k' -- the k-mer size */
        size_t ksize() const;

        /** normalize row i, so that the sum adds to 1, specifying a probability
         * distribution over the k-mers. */
        void make_distribution(size_t i);

        /** normalize every row to form a distribution */
        void make_distribution();

        /** normalize row i, so that it represents the distribution over
         * nucleotide u within the h-mer, conditioned on all the other
         * positions.
         *
         * */
        void make_conditional_distribution(size_t i, size_t u, size_t h);

        /** log-transform a row, assuming it specifies a distribution of h-mers */
        void dist_log_transform_row(size_t i, int h);

    private:
        size_t k;      //< size of k-mer
        size_t size1;  //< number of rows
        size_t size2;  //< number of columns (4^k)
        double* A;     //< data
};



#endif



