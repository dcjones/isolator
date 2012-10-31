/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#include <cmath>
#include <cstring>

#include "kmer_matrix.hpp"
#include "logger.hpp"


/* return 4^x */
static unsigned int four_pow(unsigned int x)
{
    return 1 << (2 * x);
}


#if 0
kmer_matrix::kmer_matrix(const YAML::Node& node)
{
    unsigned int k_;
    node["k"] >> k_;
    k = (size_t) k;

    unsigned int size1_;
    node["size1"] >> size1_;
    size1 = (size_t) size1_;

    unsigned int size2_;
    node["size2"] >> size2_;
    size2 = (size_t) size2_;


    A = new double[size1 * size2];

    const YAML::Node& node_A = node["A"];
    size_t i;
    for (i = 0; i < size1 * size2; ++i) {
        node_A[i] >> A[i];
    }
}
#endif


kmer_matrix::kmer_matrix(size_t size1, size_t k)
    : k(k)
    , size1(size1)
{
    size2 = (size_t) four_pow((unsigned int) k);
    A = new double[size1 * size2];
}


kmer_matrix::kmer_matrix(const kmer_matrix& other)
    : k(other.k)
    , size1(other.size1)
    , size2(other.size2)
{
    A = new double[size1 * size2];
    memcpy(A, other.A, size1 * size2 * sizeof(double));
}


kmer_matrix::~kmer_matrix()
{
    delete [] A;
}


void kmer_matrix::operator = (const kmer_matrix& other)
{

    if (size1 != other.size1 || size2 != other.size2) {
        size1 = other.size1;
        size2 = other.size2;
        k     = other.k;
        delete [] A;
        A = new double[size1 * size2];
    }

    memcpy(A, other.A, size1 * size2 * sizeof(double));
}

#if 0
void kmer_matrix::to_yaml(YAML::Emitter& out) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "k";
    out << YAML::Value << (unsigned int) k;

    out << YAML::Key   << "size1";
    out << YAML::Value << (unsigned int) size1; 

    out << YAML::Key   << "size2";
    out << YAML::Value << (unsigned int) size2; 

    out << YAML::Key   << "A";
    out << YAML::Flow;
    out << YAML::Value;
    out << YAML::BeginSeq;

    size_t i;
    for (i = 0; i < size1 * size2; ++i) {
        out << A[i];
    }

    out << YAML::EndSeq;

    out << YAML::EndMap;
}
#endif


double& kmer_matrix::operator()(size_t i, size_t j)
{
    return A[i * size2 + j];
}




void kmer_matrix::set_all(double x)
{
    const size_t N = size1 * size2;
    size_t a;
    for (a = 0; a < N; ++a) A[a] = x;
}

void kmer_matrix::set_row(size_t i, double x)
{
    double* B = A + (i * size2);
    size_t j;
    for (j = 0; j < size2; ++j) {
        B[j] = x;
    }
}

void kmer_matrix::get_row(size_t i, double* xs)
{
    double* B = A + (i * size2);
    size_t j;
    for (j = 0; j < size2; ++j) {
        xs[j] = B[j];
    }
}

void kmer_matrix::set_row(size_t i, double* xs)
{
    double* B = A + (i * size2);
    size_t j;
    for (j = 0; j < size2; ++j) B[j] = xs[j];
}


size_t kmer_matrix::nrows() const
{
    return size1;
}

size_t kmer_matrix::ncols() const
{
    return size2;
}


size_t kmer_matrix::ksize() const
{
    return k;
}


void kmer_matrix::make_distribution(size_t i)
{
    double* B = A + (i * size2);
    double z = 0.0;
    size_t j;
    for (j = 0; j < size2; ++j) z += B[j];
    for (j = 0; j < size2; ++j) B[j] /= z;
}


void kmer_matrix::make_distribution()
{
    size_t i;
    for (i = 0; i < size1; ++i) make_distribution(i);
}


void kmer_matrix::make_conditional_distribution(size_t i, size_t u, size_t h)
{
    /* row */
    double* B = A + (i * size2);

    /* low order kmer */
    kmer L;
    kmer L_max = four_pow(h - u - 1);

    /* high order kmer */
    kmer H;
    kmer H_max = four_pow(u);

    kmer K;
    kmer nt;

    double z;

    for (H = 0; H < H_max; ++H) {
        for (L = 0; L < L_max; ++L) {

            z = 0.0;
            for (nt = 0; nt < 4; ++nt) {
                K = (H << (2 * (h - u))) | (nt << (2 * (h - u - 1))) | L;
                z += B[K];
            }

            for (nt = 0; nt < 4; ++nt) {
                K = (H << (2 * (h - u))) | (nt << (2 * (h - u - 1))) | L;
                B[K] /= z;
            }
        }
    }
}

void kmer_matrix::dist_log_transform_row(size_t i, int h)
{
    double* B = A + (i * size2);
    size_t N = four_pow(h);
    size_t j;
    for (j = 0; j < N; ++j) {
        B[j] = log(B[j]);
    }
}

