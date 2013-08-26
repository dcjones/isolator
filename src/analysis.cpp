
#include <gsl/gsl_statistics_double.h>
#include <boost/foreach.hpp>
#include <vector>

#include "analysis.hpp"

using namespace std;
using namespace boost::numeric::ublas;


/* Blah. Is this going to work at all?
 *
 * There's not really a strong shrinkage effect for means across conditions, so
 * this will likely be dominated by genes with low expression, right?
 *
 * Well, let's try it, and if we need to add more levels in the goddamn
 * hiearchy, then so be it.
 */

void AvgPairwiseAbsLog2Fc::run(const std::vector<matrix<float> >& mu_samples,
                               const std::vector<set<GeneID> >& tss_gene_ids,
                               const std::vector<set<TranscriptID> >& tss_transcript_ids)
{
    if (mu_samples.empty()) return;

    FILE* out = fopen("avgabslog2fc.tsv", "w");
    if (!out) {
        Logger::abort("Can't open avgabslog2fc.tsv for writing.");
    }
    fprintf(out, "gene_ids\ttranscript_ids\tavgabslog2fc\n");

    size_t n_cond = mu_samples[0].size1();
    size_t n_tss  = mu_samples[0].size2();

    std::vector<double> abslog2fc(mu_samples.size());

    for (unsigned int j = 0; j < n_tss; ++j) {
        for (unsigned int u = 0; u < n_cond; ++u) {
            for (unsigned int v = u + 1; v < n_cond; ++v) {
                for (unsigned int k = 0; k < mu_samples.size(); ++k) {
                    abslog2fc[k] = abs(mu_samples[k](u, j) - mu_samples[k](u, j));
                }
            }
        }

        bool first = true;
        BOOST_FOREACH (const GeneID& gene_id, tss_gene_ids[j]) {
            if (!first) {
                fprintf(out, ",");
            }
            fprintf(out, "%s", gene_id.get().c_str());
            first = false;
        }
        fprintf(out, "\t");

        first = true;
        BOOST_FOREACH (const TranscriptID& transcript_id, tss_transcript_ids[j]) {
            if (!first) {
                fprintf(out, ",");
            }
            fprintf(out, "%s", transcript_id.get().c_str());
            first = false;
        }
        fprintf(out, "\t%e\n", gsl_stats_mean(&abslog2fc.at(0), 1, abslog2fc.size()));
    }

    fclose(out);
}


