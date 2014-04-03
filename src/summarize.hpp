
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

#include "hdf5.hpp"
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
    std::vector<std::string> sample_names;
    unsigned int rng_seed;
};


enum GeneFeatureType {
    GENE_FEATURE_CASSETTE_EXON,
    GENE_FEATURE_RETAINED_INTRON
};


class Summarize
{
    public:
        Summarize(const char* filename);
        ~Summarize();

        // pre-baked summarizations
        void transcript_expression(FILE* output,
                                   double credible_interval,
                                   bool unnormalized, bool splicing_rate);

        void gene_expression(FILE* output,
                             double credible_interval,
                             bool unnormalized);

        void differential_transcription(FILE* output, double credible_interval,
                                        double effect_size);

        void differential_splicing(FILE* output, double credible_interval,
                                   double effect_size);

        void differential_feature_splicing(FILE* output,
                                           double credible_interval,
                                           double effect_size);

        void condition_splicing(FILE* output, double credible_interval);

        void condition_splicing_sigma(FILE* output, double credible_interval);


        void experiment_splicing(FILE* output, double credible_interval);

        void experiment_splicing_sigma(FILE* output, double credible_interal);

        // get raw sampler output for a particular trasncript or tgroup
        std::vector<float> transcript_experiment_splicing(const char* transcript_id);
        std::vector<std::vector<float> > transcript_condition_splicing(const char* transcript_id);


        // TODO: these are not exposed. Either get rid of them, or make them
        // accessable.
        void median_condition_tgroup_expression(FILE* output);

        void median_experiment_tgroup_sd(FILE* output);

        void tgroup_fold_change(FILE* output,
                                unsigned int condition_a,
                                unsigned int condition_b);

        void expression_samples(FILE* output);
        void condition_pairwise_splicing(FILE* output);
        //void cassette_exon_pairwise_splicing(FILE* output);

        void read_metadata(IsolatorMetadata& metadata);

        const std::vector<TranscriptID>& get_transcript_ids();
        const std::vector<GeneID>& get_gene_ids();
        const std::vector<unsigned int>& get_tgroups();

    private:
        void median_ci_transcript_expression(
                boost::numeric::ublas::matrix<float>* med,
                boost::numeric::ublas::matrix<float>* lower,
                boost::numeric::ublas::matrix<float>* upper,
                double credible_interval, bool unnormalized,
                bool splicing_rate);

        void median_ci_gene_expression(
                boost::numeric::ublas::matrix<float>* med,
                boost::numeric::ublas::matrix<float>* lower,
                boost::numeric::ublas::matrix<float>* upper,
                double credible_interval, bool unnormalized);

        void read_condition_tgroup_mean(unsigned int condition,
                                        boost::numeric::ublas::matrix<float>& data);

        void condition_splicing(std::vector<boost::multi_array<float, 3> >& output);

        void read_gene_features(
                std::vector<Interval>& feature_intervals,
                std::vector<std::vector<unsigned int> >& including_tids,
                std::vector<std::vector<unsigned int> >& excluding_tids,
                std::vector<GeneFeatureType>& feature_types);

        void read_tgroup_mean(boost::multi_array<float, 3>& output);

        hid_t h5_file;

        // number of sampler iterations
        size_t num_samples;

        // number of samples
        size_t K;

        // number of conditions
        size_t C;

        // number of transcripts
        size_t N;

        // number of tgroups
        size_t T;

        IsolatorMetadata metadata;

        std::vector<std::string> condition_names;

        // transcript_id indexed by tid
        std::vector<TranscriptID> transcript_ids;

        // gene_name indexed by tid
        std::vector<GeneName> gene_names;

        // gene_id indexed by tid
        std::vector<GeneID> gene_ids;

        // tgroup indexed by tid
        std::vector<unsigned int> tgroup;

        // tids belonging to each tgroup
        std::vector<std::vector<unsigned int> > tgroup_tids;

        // indexes of tgroups with multiple transcripts
        std::vector<unsigned int> spliced_tgroup_indexes;

        // normalization factor indexed by sample_num, sample
        boost::numeric::ublas::matrix<double> scale;
};


#endif
