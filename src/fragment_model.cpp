#include <boost/foreach.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/unordered_map.hpp>
#include <cstdio>

#include "constants.hpp"
#include "fastmath.hpp"
#include "fragment_model.hpp"
#include "logger.hpp"
#include "queue.hpp"
#include "read_set.hpp"
#include "samtools/sam.h"
#include "samtools/samtools_extra.h"

extern "C" {
#include "samtools/khash.h"
KHASH_MAP_INIT_STR(s, int)
static void supress_khash_unused_function_warnings() __attribute__ ((unused));
static void supress_khash_unused_function_warnings()
{
    UNUSED(kh_init_s);
    UNUSED(kh_destroy_s);
    UNUSED(kh_put_s);
    UNUSED(kh_clear_s);
    UNUSED(kh_del_s);
}
}


/* An interval used for fragment model parameter estimation. */
class FragmentModelInterval
{
    public:
        enum IntervalType
        {
            INTERGENIC,
            CONSENSUS_EXONIC,
            EXONIC,
            THREE_PRIME
        } type;

        FragmentModelInterval(const Interval& interval,
                              IntervalType type)
            : type(type)
            , seqname(interval.seqname)
            , start(interval.start)
            , end(interval.end)
            , strand(interval.strand)
            , tpdist(-1)
            , tlen(-1)
            , tid(-1)
        {
        }

        FragmentModelInterval(SeqName seqname,
                              pos_t start,
                              pos_t end,
                              strand_t strand,
                              IntervalType type)
            : type(type)
            , seqname(seqname)
            , start(start)
            , end(end)
            , strand(strand)
            , tpdist(-1)
            , tlen(-1)
            , tid(-1)
        {
        }

        void add_alignment(long idx, const bam1_t* b)
        {
            rs.add_alignment(idx, b);
        }

        bool operator < (const FragmentModelInterval& other) const
        {
            if      (tid != other.tid)     return tid < other.tid;
            else if (start != other.start) return start < other.start;
            else                           return end < other.end;
        }

        void clear()
        {
            rs.clear();
        }

        ReadSet rs;

        SeqName seqname;
        pos_t start, end;
        strand_t strand;

        // distance this interval is from the 3' end of the transcript,
        // if applicable.
        pos_t tpdist;

        // transcript lengt,h if applicable
        pos_t tlen;

        /* A the sequence ID assigned by the BAM file, so we can arrange
         * intervals in the same order. */
        int32_t tid;
};


struct FragmentModelIntervalPtrCmp
{
    bool operator () (const FragmentModelInterval* a,
                      const FragmentModelInterval* b)
    {
        return *a < *b;
    }
};


void sam_scan(std::vector<FragmentModelInterval*>& intervals,
              AlnIndex& alnindex,
              Queue<FragmentModelInterval*>& q,
              std::set<std::string> excluded_seqs,
              std::vector<FragmentModelInterval*>& seqbias_intervals,
              PosTable& seqbias_sense_pos,
              PosTable& seqbias_antisense_pos,
              PosTable& gcbias_pos,
              const char* bam_fn,
              const char* task_name)
{
    /* Measure file size to monitor progress. */
    size_t input_size = 0;
    size_t input_block_size = 1000000;
    if (task_name) {
        FILE* f = fopen(bam_fn, "rb");
        if (f) {
            fseek(f, 0, SEEK_END);
            input_size = (size_t) ftell(f);
            fclose(f);
        }
        Logger::push_task(task_name, input_size / input_block_size);
    }

    samfile_t* bam_f;
    bam_f = samopen(bam_fn, "rb", NULL);
    if (bam_f == NULL) {
        bam_f = samopen(bam_fn, "r", NULL);
    }
    if (bam_f == NULL) {
        Logger::abort("Can't open SAM/BAM file %s.\n", bam_fn);
    }

    /* Sort the intervals in the same order as the BAM file. */
    bam_init_header_hash(bam_f->header);
    khash_t(s)* tbl = reinterpret_cast<khash_t(s)*>(bam_f->header->hash);

    std::vector<FragmentModelInterval*>::iterator i;
    khiter_t k;
    for (i = intervals.begin(); i != intervals.end(); ++i) {
        k = kh_get(s, tbl, (*i)->seqname.get().c_str());
        if (k == kh_end(tbl)) (*i)->tid = -1;
        else (*i)->tid = kh_value(tbl, k);
    }
    std::sort(intervals.begin(), intervals.end(), FragmentModelIntervalPtrCmp());

    for (i = seqbias_intervals.begin(); i != seqbias_intervals.end(); ++i) {
        k = kh_get(s, tbl, (*i)->seqname.get().c_str());
        if (k == kh_end(tbl)) (*i)->tid = -1;
        else (*i)->tid = kh_value(tbl, k);
    }
    std::sort(seqbias_intervals.begin(), seqbias_intervals.end(), FragmentModelIntervalPtrCmp());

    /* First interval which the current read may be contained in. */
    size_t j, j0 = 0;
    size_t n = intervals.size();

    /* First seqbias interval in which the current read may be contained. */
    std::vector<FragmentModelInterval*>::iterator sbi = seqbias_intervals.begin();

    size_t last_file_pos = 0, file_pos;
    size_t read_num = 0;
    bool excluded_seq = false;

    /* Read the reads. */
    bam1_t* b = bam_init1();
    int32_t last_tid = -1;
    int32_t last_pos = -1;
    while (samread(bam_f, b) >= 0) {
        ++read_num;

        if (read_num % 1000 == 0) {
            file_pos = samtell(bam_f);
            if (file_pos >= last_file_pos + input_block_size && input_size > 0) {
                Logger::get_task(task_name).inc();
                last_file_pos = file_pos;
            }
        }

        if (b->core.flag & BAM_FUNMAP || b->core.tid < 0) continue;

        if (b->core.tid != last_tid) {
            std::string seqname(bam_f->header->target_name[b->core.tid]);
            excluded_seq = excluded_seqs.find(seqname) != excluded_seqs.end();
        }

        if (b->core.tid < last_tid ||
            (b->core.tid == last_tid && b->core.pos < last_pos)) {
            Logger::abort(
                    "The input SAM/BAM file must be sorted. "
                    "Please run: 'samtools sort'.");
        }

        last_tid = b->core.tid;
        last_pos = b->core.pos;

        if (excluded_seq) continue;

        long idx = alnindex.add(bam1_qname(b));

        // Process seqbias intervals
        while (sbi != seqbias_intervals.end() && (*sbi)->tid < b->core.tid) ++sbi;
        while (sbi != seqbias_intervals.end() && (*sbi)->tid == b->core.tid &&
               (*sbi)->end < b->core.pos) ++sbi;

        // Check the next two intervals, as it might overlap one interval on
        // both strands.
        for (std::vector<FragmentModelInterval*>::iterator sbj = sbi;
             sbj != seqbias_intervals.end() && sbj != sbi + 2; ++sbj) {
            pos_t start = (pos_t) bam_truepos(&b->core, bam1_cigar(b)),
                  end = (pos_t) bam_trueend(&b->core, bam1_cigar(b)) - 1;
            if ((*sbj)->tid == b->core.tid &&
                (*sbj)->start <= start && end <= (*sbj)->end) {
                strand_t b_strand = bam1_strand(b) ? strand_neg : strand_pos;
                pos_t pos = b_strand ? end : start;
                if ((*sbj)->strand == b_strand) {
                    seqbias_sense_pos.add(b->core.tid, pos, bam1_strand(b),
                                          (*sbj)->start, (*sbj)->end, bam_f);
                }
                else {
                    seqbias_antisense_pos.add(b->core.tid, pos, bam1_strand(b),
                                              (*sbj)->start, (*sbj)->end, bam_f);
                }

                if ((*sbj)->end - (*sbj)->start + 1 >= constants::gcbias_min_seq_len) {
                    gcbias_pos.add(b->core.tid, pos, bam1_strand(b),
                                   (*sbj)->start, (*sbj)->end, bam_f);
                }
            }
        }

        // Add reads to intervals in which they are contained.
        for (j = j0; j < n; ++j) {
            if (b->core.tid < intervals[j]->tid) break;
            if (b->core.tid > intervals[j]->tid) {
                assert(j == j0);
                q.push(intervals[j0++]);
                continue;
            }

            pos_t pos = bam_truepos(&b->core, bam1_cigar(b));
            if (pos < intervals[j]->start) break;
            if (pos > intervals[j]->end) {
                if (j == j0) {
                    q.push(intervals[j0++]);
                }
                continue;
            }

            pos_t b_end = (pos_t) bam_trueend(&b->core, bam1_cigar(b)) - 1;
            if (b_end <= intervals[j]->end) {
                intervals[j]->add_alignment(idx, b);
            }
        }
    }

    for(; j0 < n; ++j0) q.push(intervals[j0]);

    bam_destroy1(b);
    samclose(bam_f);

    Logger::debug("alignement index is %0.2fMB",
                  (double) alnindex.used_memory() / 1e6);

    if (task_name) Logger::pop_task(task_name);
}


/* A thread used to accumulate statistics to estimate the parameters of
 * the fragment model. */
class FragmentModelThread
{
    public:
        FragmentModelThread(Queue<FragmentModelInterval*>& q)
            : q(q)
            , t(NULL)
        {
            strand_bias[0] = strand_bias[1] = 0;
        }

        ~FragmentModelThread()
        {
            delete t;
        }

        void run()
        {
            FragmentModelInterval* interval;
            ReadSet::UniqueReadCounts read_counts;

            while (true) {
                interval = q.pop();
                if (interval == NULL) break;

                /* dispatch! */
                if (interval->type == FragmentModelInterval::INTERGENIC) {
                    /* TODO: measure additive noise */
                }
                else if (interval->type == FragmentModelInterval::CONSENSUS_EXONIC) {
                    interval->rs.make_unique_read_counts(read_counts);
                    measure_fragment_lengths(read_counts);
                    measure_strand_bias(interval->strand, read_counts);
                    read_counts.clear();
                }
                else if (interval->type == FragmentModelInterval::EXONIC) {
                    interval->rs.make_unique_read_counts(read_counts);
                    measure_tpbias(interval->start, interval->end,
                                   interval->strand, interval->tpdist,
                                   interval->tlen, read_counts);
                    read_counts.clear();
                }

                delete interval;
            }
        }

        void start()
        {
            if (t != NULL) return;
            t = new boost::thread(boost::bind(&FragmentModelThread::run, this));
        }

        void join()
        {
            t->join();
            delete t;
            t = NULL;
        }

        /* strand_bias[0] counts the number of times the read's strand agrees
         * with the transcripts, and strand_bias[1] the number of disagreements.
         * */
        unsigned long strand_bias[2];

        /* A hash table mapping fragment lengths to number of observations */
        boost::unordered_map<unsigned int, unsigned int> frag_lens;

        /* A hash table mapping a number k to the number of intergenic
         * positions with k read starts. */
        boost::unordered_map<unsigned int, unsigned int> noise_counts;

        /* Mate1 sequence bias across the current interval. */
        std::vector<double> mate1_seqbias[2];

        /* (tpdist, tlen) pairs used to fit three prime bias model. */
        std::vector<std::pair<pos_t, pos_t> > tpbias_pairs;

    private:
        void measure_fragment_lengths(const ReadSet::UniqueReadCounts& counts)
        {
            boost::unordered_map<unsigned int, unsigned int>::iterator c;
            ReadSet::UniqueReadCounts::const_iterator i;
            for (i = counts.begin(); i != counts.end(); ++i) {
                AlignedReadIterator j(*i->first);
                for (; j != AlignedReadIterator(); ++j) {
                    if (j->mate1 == NULL || j->mate2 == NULL) continue;
                    if (j->mate1->cigar_len != 1 || j->mate2->cigar_len != 1) continue;

                    pos_t len = j->naive_frag_len();

                    if (len <= 0 || len > constants::max_frag_len) continue;

                    c = frag_lens.find(len);
                    if (c == frag_lens.end()) {
                        frag_lens.insert(std::make_pair((unsigned int) len, 1));
                    }
                    else {
                        c->second += 1;
                    }
                }
            }
        }

        void measure_strand_bias(strand_t strand,
                                 const ReadSet::UniqueReadCounts& counts)
        {
            ReadSet::UniqueReadCounts::const_iterator i;
            for (i = counts.begin(); i != counts.end(); ++i) {
                AlignedReadIterator j(*i->first);
                for (; j != AlignedReadIterator(); ++j) {
                    if (!j->mate1) continue;
                    strand_bias[j->mate1->strand == strand ? 0 : 1]++;
                }
            }
        }

        void measure_tpbias(pos_t start, pos_t end, strand_t strand,
                            pos_t tpdist, pos_t tlen,
                            const ReadSet::UniqueReadCounts& counts)
        {
            ReadSet::UniqueReadCounts::const_iterator i;
            for (i = counts.begin(); i != counts.end(); ++i) {
                AlignedReadIterator j(*i->first);
                for (; j != AlignedReadIterator(); ++j) {
                    if (strand == strand_pos) {
                        if (j->mate1 && j->mate1->strand == 0 &&
                                start <= j->mate1->start && j->mate1->start <= end) {
                            tpbias_pairs.push_back(
                                    std::make_pair(end - j->mate1->start + tpdist, tlen));
                        }
                        else if (j->mate2 && j->mate2->strand == 0 &&
                                    start <= j->mate2->start && j->mate2->start <= end) {
                            tpbias_pairs.push_back(
                                    std::make_pair(end - j->mate2->start + tpdist, tlen));
                        }
                    }
                    else {
                        if (j->mate1 && j->mate1->strand == 1 &&
                                start <= j->mate1->end && j->mate1->end <= end) {
                            tpbias_pairs.push_back(
                                std::make_pair(j->mate1->end - start + tpdist, tlen));
                        }
                        else if (j->mate2 && j->mate2->strand == 1 &&
                                    start <= j->mate2->end && j->mate2->end <= end) {
                            tpbias_pairs.push_back(
                                std::make_pair(j->mate2->end - start + tpdist, tlen));
                        }
                    }
                }
            }
        }

        Queue<FragmentModelInterval*>& q;
        boost::thread* t;
};


FragmentModel::FragmentModel()
{
    sb[0] = sb[1] = NULL;
    frag_len_dist = NULL;
}


FragmentModel::~FragmentModel()
{
    delete sb[0];
    delete sb[1];
    delete gcbias;
    delete tpbias;
    delete fragbias;
    delete frag_len_dist;
}


void FragmentModel::estimate(TranscriptSet& ts,
        const char* bam_fn,
        const char* fa_fn,
        bool use_seqbias_correction,
        bool use_gc_correction,
        bool use_3p_correction,
        bool use_frag_correction,
        bool tabulate_bias,
        std::set<std::string> excluded_seqs,
        std::set<std::string> bias_training_seqnames)
{
    Queue<FragmentModelInterval*> q(constants::max_estimate_queue_size);

    std::vector<FragmentModelInterval*> intervals;

    // consensus exonic intervals for training fragment length distribution
    std::vector<Interval> exonic;
    ts.get_consensus_exonic(exonic);
    std::vector<Interval>::iterator interval;
    for (interval = exonic.begin(); interval != exonic.end(); ++interval) {
        intervals.push_back(new FragmentModelInterval(
                    *interval,
                    FragmentModelInterval::CONSENSUS_EXONIC));
    }
    Logger::info("%lu consensus exons", (unsigned long) intervals.size());

    // exonic intervals for training tpbias
    if (use_3p_correction) {
        size_t num_tpbias_transcripts = 0;
        for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
            pos_t offset = 0;
            pos_t tlen = t->exonic_length();

            if (tlen < constants::tpbias_min_tlen ||
                tlen > constants::tpbias_max_tlen) continue;
            if (num_tpbias_transcripts++ > constants::tpbias_max_transcripts) break;

            BOOST_FOREACH (const Exon& exon, *t) {
                intervals.push_back(
                        new FragmentModelInterval(t->seqname, exon.start, exon.end,
                                                  t->strand, FragmentModelInterval::EXONIC));
                intervals.back()->tlen = t->exonic_length();
                intervals.back()->tpdist =
                    t->strand == strand_pos ?
                        tlen - offset - (exon.end - exon.start + 1): offset;
                offset += exon.end - exon.start + 1;
            }
        }
    }

    // exonic intervals for training seqbias
    std::vector<FragmentModelInterval*> seqbias_intervals;
    if (fa_fn) {
        for (interval = exonic.begin(); interval != exonic.end(); ++interval) {
            if (!bias_training_seqnames.empty() &&
                bias_training_seqnames.find(interval->seqname.get()) == bias_training_seqnames.end()) {
                continue;
            }
            seqbias_intervals.push_back(new FragmentModelInterval(
                *interval,
                FragmentModelInterval::EXONIC));
        }
    }

    std::vector<FragmentModelThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new FragmentModelThread(q));
        threads.back()->start();
    }

    PosTable* seqbias_sense_pos = new PosTable();
    PosTable* seqbias_antisense_pos = new PosTable();
    PosTable* gcbias_pos = new PosTable();

    std::string task_name = std::string("Indexing reads (") +
                            std::string(bam_fn) +
                            std::string(")");

    sam_scan(intervals, alnindex, q, excluded_seqs,
             seqbias_intervals, *seqbias_sense_pos, *seqbias_antisense_pos,
             *gcbias_pos, bam_fn, task_name.c_str());

    // delete seqbias intervals
    for (std::vector<FragmentModelInterval*>::iterator i = seqbias_intervals.begin();
            i != seqbias_intervals.end(); ++i) {
        delete *i;
    }
    seqbias_intervals.clear();

    // train seqbias
    if (fa_fn && use_seqbias_correction) {
        std::string task_name0 = std::string("Training sense sequence bias (") +
                                 std::string(bam_fn) +
                                 std::string(")");

        std::string task_name1 = std::string("Training antisense sequence bias (") +
                                 std::string(bam_fn) +
                                 std::string(")");

        sb_tabulation[0].order = 1;
        sb[0] = new sequencing_bias(fa_fn, *seqbias_sense_pos,
                                    constants::seqbias_num_reads,
                                    constants::seqbias_left_pos,
                                    constants::seqbias_right_pos,
                                    task_name0.c_str(),
                                    tabulate_bias ? &sb_tabulation[0] : NULL);

        sb_tabulation[1].order = 1;
        sb[1] = new sequencing_bias(fa_fn, *seqbias_antisense_pos,
                                    constants::seqbias_num_reads,
                                    constants::seqbias_left_pos,
                                    constants::seqbias_right_pos,
                                    task_name1.c_str(),
                                    tabulate_bias ? &sb_tabulation[1] : NULL);
    }
    else sb[0] = sb[1] = NULL;

    delete seqbias_sense_pos;
    delete seqbias_antisense_pos;

    Logger::debug("Indexed %lu reads", alnindex.size());

    for (size_t i = 0; i < constants::num_threads; ++i) {
        q.push(NULL);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads[i]->join();
    }

    /* Collect strand specificity statistics. */
    unsigned long strand_bias[2] = {0, 0};
    for (size_t i = 0; i < constants::num_threads; ++i) {
        strand_bias[0] += threads[i]->strand_bias[0];
        strand_bias[1] += threads[i]->strand_bias[1];
    }

    strand_specificity = (float) strand_bias[0];
    strand_specificity /= (float) (strand_bias[0] + strand_bias[1]);

    double negentropy =
        strand_specificity * log2(strand_specificity) +
        (1.0 - strand_specificity) * log2(1.0 - strand_specificity);
    Logger::debug("Strand specificity: %0.1f%%", 100.0 * (1.0 + negentropy));

    /* Collect fragment length statistics. */
    std::map<unsigned int, unsigned int> frag_lens;
    boost::unordered_map<unsigned int, unsigned int>::iterator f;
    std::map<unsigned int, unsigned int>::iterator g;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        for (f = threads[i]->frag_lens.begin();
             f != threads[i]->frag_lens.end();
             ++f) {
            g = frag_lens.find(f->first);
            if (g == frag_lens.end()) {
                frag_lens.insert(*f);
            }
            else {
                g->second += f->second;
            }
        }
    }

    size_t frag_len_pe_reads = 0;
    for (g = frag_lens.begin(); g != frag_lens.end(); ++g) {
        frag_len_pe_reads += g->second;
    }

    if (frag_len_pe_reads > constants::frag_len_min_pe_reads) {
        unsigned int* frag_len_vals = new unsigned int [frag_lens.size()];
        unsigned int* frag_len_lens = new unsigned int [frag_lens.size()];
        unsigned int sum_lens = 0;

        size_t i;
        for (g = frag_lens.begin(), i = 0; g != frag_lens.end(); ++g, ++i) {
            frag_len_vals[i] = g->first;
            frag_len_lens[i] = g->second;
            sum_lens += g->second;
        }

        frag_len_dist = new EmpDist(frag_len_vals, frag_len_lens, frag_lens.size());
        {
            FILE* out = fopen("frag-len-dist.txt", "w");
            for (pos_t fl = 1; fl < 1000; ++fl) {
                fprintf(out, "%ld\t%e\t%e\n",
                        fl, (double) frag_len_dist->pdf(fl),
                        (double) frag_len_dist->cdf(fl));
            }

            fclose(out);
        }

        delete [] frag_len_vals;
        delete [] frag_len_lens;

        Logger::info("Median fragment-length: %d", (int) frag_len_dist->median());
    }
    else {
        Logger::info("Too few paired-end reads to estimate fragment length distribution.");
    }

    // train GC bias
    if (use_gc_correction && fa_fn) {
        std::string task_name = std::string("Training GC bias (") +
                                std::string(bam_fn) +
                                std::string(")");

        gcbias = new GCBias(fa_fn, *gcbias_pos, frag_len_med(), sb,
                            task_name.c_str());
    }
    else gcbias = NULL;

    delete gcbias_pos;

    // train 3' bias
    if (use_3p_correction) {
        std::vector<std::pair<pos_t, pos_t> > tpdists;
        for (size_t i = 0; i < constants::num_threads; ++i) {
            tpdists.insert(tpdists.end(),
                           threads[i]->tpbias_pairs.begin(),
                           threads[i]->tpbias_pairs.end());
        }

        tpbias = new TPBias(tpdists);
        Logger::info("3' bias: %0.3e", tpbias->p);
    }
    else tpbias = NULL;

    for (size_t i = 0; i < constants::num_threads; ++i) {
        delete threads[i];
    }

    // train fragmentation bias
    if (use_frag_correction) {
        fragbias = new FragBias(this);
    }
    else fragbias = NULL;
}


float FragmentModel::frag_len_p(pos_t frag_len)
{
    if (frag_len_dist) return frag_len_dist->pdf(frag_len);
    else return boost::math::pdf(
            boost::math::normal_distribution<double>(constants::frag_len_mu,
                                                    constants::frag_len_sd),
            frag_len);
}


float FragmentModel::frag_len_c(pos_t frag_len)
{
    if (frag_len_dist) return frag_len_dist->cdf(frag_len);
    else return boost::math::cdf(
            boost::math::normal_distribution<double>(constants::frag_len_mu,
                                                    constants::frag_len_sd),
            frag_len);
}


float FragmentModel::frag_len_med()
{
    if (frag_len_dist) return frag_len_dist->median();
    else return constants::frag_len_mu;
}


