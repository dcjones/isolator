
#include <boost/foreach.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <cstdio>
#include <algorithm>

#include "constants.hpp"
#include "gcbias.hpp"
#include "logger.hpp"
#include "seqbias/twobitseq.hpp"


static size_t gc_count(char* seq, size_t n)
{
	size_t count = 0;
	for (size_t i = 0; i < n; ++i) {
		if (seq[i] == 'g' || seq[i] == 'G' || seq[i] == 'c' || seq[i] == 'C') ++count;
	}
	return count;
}


struct ReadPosSeqnameCmp
{
    bool operator () (const ReadPos& a, const ReadPos& b)
    {
        return a.seqname < b.seqname;
    }
};


GCBias::GCBias(const char* ref_filename, PosTable& foreground_position_table,
	           pos_t median_frag_len,
               sequencing_bias* seqbias[2],
               const char* task_name)
{
	faidx_t* ref_file = fai_load(ref_filename);
	if (!ref_file) {
        Logger::abort("Can't open fasta file '%s'.", ref_filename);
	}

	std::vector<ReadPos> foreground_positions;
    const size_t max_dump = 10000000;
	foreground_position_table.dump(foreground_positions, max_dump);
	std::sort(foreground_positions.begin(), foreground_positions.end(), ReadPosSeqnameCmp());

	Logger::push_task(task_name, foreground_positions.size());
	LoggerTask& task = Logger::get_task(task_name);

    typedef std::pair<float, float> WeightedGC;

	std::vector<WeightedGC> foreground_gc, background_gc;

    int  seqlen = 0;
    SeqName curr_seqname;
    char* seq = NULL;
    twobitseq tbseq;
    twobitseq tbseqrc;
    rng_t rng;

    pos_t L = seqbias[0]->getL();

    std::vector<ReadPos>::iterator i;
    for (i = foreground_positions.begin(); i != foreground_positions.end(); ++i) {
        if (i->seqname != curr_seqname) {
            if (seq) free(seq);

            seq = faidx_fetch_seq(ref_file, i->seqname.get().c_str(), 0, INT_MAX, &seqlen);
            Logger::debug("read sequence %s.", i->seqname.get().c_str());

            if (seq == NULL) {
                Logger::warn("warning: reference sequence not found, skipping.");
            }
            else {
                for (char* c = seq; *c; c++) *c = tolower(*c);
                tbseq = seq;
                tbseqrc = tbseq;
                tbseqrc.revcomp();
            }

            curr_seqname = i->seqname;
        }

        if (seq == NULL) continue;

        if (i->count > 16 || i->count < 3) continue;

        // sample background position
        boost::random::uniform_int_distribution<pos_t> random_uniform(
                i->start, i->end - median_frag_len + 1);
        pos_t pos = random_uniform(rng);
        float gc = (float) gc_count(seq + pos, median_frag_len) / median_frag_len;
        float sb = seqbias[0]->get_bias(tbseq, pos - L) *
                   seqbias[1]->get_bias(tbseqrc, seqlen - pos - 1 - L);
        background_gc.push_back(WeightedGC(gc, 1.0 / sb));

        // sample foreground position
        if (i->strand == 0) {
            if (i->pos >= i->start && i->pos + median_frag_len - 1 <= i->end) {
                float sb = seqbias[0]->get_bias(tbseq, i->pos - L) *
                           seqbias[1]->get_bias(tbseqrc, seqlen - i->pos - 1 - L);
                foreground_gc.push_back(
                    WeightedGC((float) gc_count(seq + i->pos, median_frag_len) / median_frag_len,
                               1.0 / sb));
            }
        } else {
            if (i->pos - median_frag_len >= i->start && i->pos <= i->end) {
                float sb = seqbias[0]->get_bias(tbseq, i->pos - median_frag_len - L) *
                           seqbias[1]->get_bias(tbseqrc, seqlen - i->pos - median_frag_len - 1 - L);
                foreground_gc.push_back(
                    WeightedGC((float) gc_count(seq + i->pos - median_frag_len, median_frag_len) / median_frag_len,
                               1.0 /sb));
            }
        }

        task.inc();
    }

#if 0
    FILE* out = fopen("gcbias.tsv", "w");
    fprintf(out, "group\tgc\n");
    BOOST_FOREACH (float value, foreground_gc) {
    	fprintf(out, "foreground\t%f\n", value);
    }
    BOOST_FOREACH (float value, background_gc) {
    	fprintf(out, "background\t%f\n", value);
    }
    fclose(out);
#endif

    // bin into equally sized bins, measure bias

    std::sort(foreground_gc.begin(), foreground_gc.end());
    std::sort(background_gc.begin(), background_gc.end());

    bins.resize(constants::gcbias_num_bins);
    for (size_t i = 0; i < constants::gcbias_num_bins; ++i) {
    	bins[i] = foreground_gc[(i+1) * (foreground_gc.size() / constants::gcbias_num_bins)].first;
    }
    bins.back() = 1.0;

    double fore_total_weight = 0.0, back_total_weight = 0.0;
    BOOST_FOREACH (const WeightedGC& gc, foreground_gc) {
        fore_total_weight += gc.second;
    }

    BOOST_FOREACH (const WeightedGC& gc, background_gc) {
        back_total_weight += gc.second;
    }

    bin_bias.resize(constants::gcbias_num_bins);
    size_t j_fore = 0, j_back = 0;
    for (size_t i = 0; i < constants::gcbias_num_bins; ++i) {
    	double fore_weight = 0.0;
    	while (j_fore < foreground_gc.size() && foreground_gc[j_fore].first < bins[i]) {
    		++j_fore;
            fore_weight += foreground_gc[j_fore].second;
    	}
        fore_weight = fore_weight / fore_total_weight;

        double back_weight = 0.0;
    	while (j_back < background_gc.size() && background_gc[j_back].first < bins[i]) {
            ++j_back;
            back_weight += background_gc[j_back].second;
        }
        back_weight = back_weight / back_total_weight;

    	bin_bias[i] = fore_weight / back_weight;

        if (!boost::math::isfinite(bin_bias[i])) bin_bias[i] = 1.0;
        if (bin_bias[i] > constants::gcbias_max_bias) {
            bin_bias[i] = constants::gcbias_max_bias;
        }
        if (bin_bias[i] < 1.0 / constants::gcbias_max_bias) {
            bin_bias[i] = 1.0 / constants::gcbias_max_bias;
        }
    }

    Logger::pop_task(task_name);
}


double GCBias::get_bias(double gc)
{
    std::vector<double>::iterator i =
        std::lower_bound(bins.begin(), bins.end(), gc);

    if (i == bins.end()) return bin_bias.back();

    return bin_bias[i - bins.begin()];
}


GCBias::~GCBias()
{

}
