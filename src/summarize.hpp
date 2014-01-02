
#ifndef ISOLATOR_SUMMARIZE_HPP
#define ISOLATOR_SUMMARIZE_HPP

// Read and produce plain-text tables from the HDF5 output generated
// by "analyze".

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "intervals.hpp"


// Data that 'isolator analyze' writes as an attribute to the hdf5 file.
struct IsolatorMetadata
{
    std::string command_line;
    std::string version;
    std::string date;
    std::string elapsed_seconds;
    std::vector<std::string> sample_filenames;
    std::vector<std::string> sample_conditions;
};


class Summarize
{
    public:
        Summarize(const char* filename);
        ~Summarize();

        // pre-baked summarizations
        void median_condition_tgroup_expression(FILE* output);
        void median_experiment_tgroup_sd(FILE* output);

        void median_transcript_expression(FILE* output);
        void median_gene_expression(FILE* output);

        void tgroup_fold_change(FILE* output,
                                unsigned int condition_a,
                                unsigned int condition_b);

        void condition_splicing(FILE* output);
        void expression_samples(FILE* output);
        void condition_pairwise_splicing(FILE* output);
        void cassette_exon_pairwise_splicing(FILE* output);

        void read_metadata(IsolatorMetadata& metadata);

    private:
        void median_transcript_expression(
                boost::numeric::ublas::matrix<float>& Q);

        void condition_splicing(std::vector<boost::multi_array<float, 3> >& output);

        void read_cassette_exons(
                std::vector<Interval>& cassette_exons,
                std::vector<std::vector<unsigned int> >& including_tids,
                std::vector<std::vector<unsigned int> >& excluding_tids);

        void read_tgroup_mean(boost::multi_array<float, 3>& output);

        hid_t h5_file;

        // number of samples
        size_t K;

        // number of transcripts
        size_t N;

        // number of tgroups
        size_t T;

        // transcript_id indexed by tid
        std::vector<std::string> transcript_ids;

        // gene_id indexed by tid
        std::vector<std::string> gene_ids;

        // tgroup indexed by tid
        std::vector<unsigned int> tgroup;

        // tids belonging to each tgroup
        std::vector<std::vector<unsigned int> > tgroup_tids;

        // indexes of tgroups with multiple transcripts
        std::vector<unsigned int> spliced_tgroup_indexes;
};


#endif
