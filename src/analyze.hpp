
#ifndef ISOLATOR_ANALYZE_HPP
#define ISOLATOR_ANALYZE_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <set>
#include <vector>

#include "sample_db.hpp"
#include "transcripts.hpp"


class Analyze
{
    public:
		Analyze(TranscriptSet& ts);
		~Analyze();

		// Add a replicate under a particular condition
		void add_sample(const char* condition_name, SampleDB& sample_data);

        void run();

    private:
        void load_quantification_data();
        void compute_depth();
        void choose_kde_bandwidth();

		TranscriptSet& ts;

		// map a condition name to data from some number of replicates
		typedef std::map<std::string, std::vector<SampleDB*> > CondMap;
		typedef std::pair<const std::string, std::vector<SampleDB*> > CondMapItem;
		CondMap data;

        // tid to transcript_ids and gene_ids
        std::vector<TranscriptID> transcript_ids;
        std::vector<GeneID> gene_ids;

        // organize by tss
		typedef std::map<Interval, std::vector<unsigned int> > TSMap;
		typedef std::pair<const Interval, std::vector<unsigned int> > TSMapItem;
		TSMap tss_group;
		TSMap tss_tts_group;

        // Transcript and gene IDs indexed by tss index.
        std::vector<std::set<GeneID> > tss_gene_ids;
        std::vector<std::set<TranscriptID> > tss_transcript_ids;
        std::vector<std::vector<unsigned int> > tss_tids;

		// map each transcript index to it's tss index
		std::vector<int> tss_index;

        // Matrix containing sampled transcript abundances from
        // 'isolator quantify'
		// Indexed by: replicate -> transcript -> sample_num
        std::vector<boost::numeric::ublas::matrix<float> > quantification;

		// condition index to sample indexes
		std::vector<std::vector<unsigned int> > condition_samples;

        // bandwidth for the kde estimator of each sample and transcript
        boost::numeric::ublas::matrix<float> bandwidth;

        // normalization constant for each sample
        std::vector<double> depth;

        // number of sequenced samples
        unsigned int K;

        // number of conditions
        unsigned int C;

        // number of transcripts
        unsigned int N;

        // number of samples generated in the quantification phase
        unsigned int M;

        // number of transcription groups (TSS groups, typicall)
        unsigned int T;

        // used to load data into the sampler
        friend class AnalyzeSamplerData;
};


#endif

