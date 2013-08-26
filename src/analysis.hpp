
#ifndef ISOLATOR_ANALYSIS_HPP
#define ISOLATOR_ANALYSIS_HPP

#include "common.hpp"
#include <vector>
#include <set>
#include <boost/numeric/ublas/matrix.hpp>

class Analysis
{
    public:
        // TODO: What does a sample look like?
        virtual void run(const std::vector<boost::numeric::ublas::matrix<float> >& mu_samples,
                         const std::vector<std::set<GeneID> >& tss_gene_ids,
                         const std::vector<std::set<TranscriptID> >& tss_transcript_ids) = 0;
};


class AvgPairwiseAbsLog2Fc : public Analysis
{
    public:
        void run(const std::vector<boost::numeric::ublas::matrix<float> >& mu_samples,
                 const std::vector<std::set<GeneID> >& tss_gene_ids,
                 const std::vector<std::set<TranscriptID> >& tss_transcript_ids);
};



#endif

