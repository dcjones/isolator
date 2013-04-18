
#include <boost/unordered_map.hpp>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_randist.h>

#include "constants.hpp"
#include "fragment_model.hpp"
#include "logger.hpp"
#include "queue.hpp"
#include "read_set.hpp"
#include "samtools/sam.h"
#include "samtools/samtools_extra.h"

extern "C" {
#include "samtools/khash.h"
KHASH_MAP_INIT_STR(s, int)
}


AlnCountTrie::AlnCountTrie()
{
    t = hattrie_create();
}


AlnCountTrie::~AlnCountTrie()
{
    hattrie_free(t);
}


void AlnCountTrie::inc_mate1(const char* id)
{
    unsigned long* val = hattrie_get(t, id, strlen(id));
    unsigned long cnt = ((*val & 0xffff) + 1) & 0xffff;
    *val = (*val & 0xffff0000) | cnt;
}


void AlnCountTrie::inc_mate2(const char* id)
{
    unsigned long* val = hattrie_get(t, id, strlen(id));
    unsigned long cnt = ((*val >> 16) + 1) & 0xffff;
    *val = (cnt << 16) | (*val & 0xffff);
}


MateCount AlnCountTrie::get(const char* id) const
{
    unsigned long* val = hattrie_tryget(t, id, strlen(id));
    return val == NULL ? MateCount(0, 0) :
                         MateCount(*val & 0xffff, (*val >> 16) & 0xffff);
}


void AlnCountTrie::set(const char* id, const MateCount& count)
{
    unsigned long* val = hattrie_get(t, id, strlen(id));
    *val = (count.first & 0xffff) | ((count.second & 0xffff) << 16);
}


bool AlnCountTrie::has(const char* id) const
{
    return hattrie_tryget(t, id, strlen(id)) != NULL;
}


size_t AlnCountTrie::size() const
{
    return hattrie_size(t);
}


AlnCountTrieIterator::AlnCountTrieIterator()
    : it(NULL)
{

}


AlnCountTrieIterator::AlnCountTrieIterator(const AlnCountTrie& t)
{
    it = hattrie_iter_begin(t.t, false);
    if (!hattrie_iter_finished(it)) {
        x.first = hattrie_iter_key(it, NULL);
        unsigned long* val = hattrie_iter_val(it);
        x.second.first  = *val & 0xffff;
        x.second.second = (*val >> 16) & 0xffff;
    }
}


AlnCountTrieIterator::~AlnCountTrieIterator()
{
    hattrie_iter_free(it);
}


void AlnCountTrieIterator::increment()
{
    hattrie_iter_next(it);
    if (!hattrie_iter_finished(it)) {
        x.first = hattrie_iter_key(it, NULL);
        unsigned long* val = hattrie_iter_val(it);
        x.second.first  = *val & 0xffff;
        x.second.second = (*val >> 16) & 0xffff;
    }
}


bool AlnCountTrieIterator::equal(const AlnCountTrieIterator& other) const
{
    if (it == NULL || hattrie_iter_finished(it)) {
        return other.it == NULL || hattrie_iter_finished(other.it);
    }
    else if (other.it == NULL || hattrie_iter_finished(other.it)) {
        return false;
    }
    else return hattrie_iter_equal(it, other.it);
}


const std::pair<const char*, MateCount>& AlnCountTrieIterator::dereference() const
{
    return x;
}


AlnIndex::AlnIndex()
{
    t = hattrie_create();
}


AlnIndex::~AlnIndex()
{
    hattrie_free(t);
}


size_t AlnIndex::size() const
{
    return hattrie_size(t);
}


void AlnIndex::clear()
{
    hattrie_clear(t);
}


void AlnIndex::add(const char* key)
{
    boost::lock_guard<boost::mutex> lock(mut);

    value_t* val = hattrie_get(t, key, strlen(key));
    if (*val == 0) {
        *val = 1 + hattrie_size(t);
    }
}


int AlnIndex::get(const char* key)
{
    boost::lock_guard<boost::mutex> lock(mut);

    value_t* val = hattrie_tryget(t, key, strlen(key));
    if (val == NULL) return -1;
    else {
        return *val - 1;;
    }
}


/* An interval used for fragment model parameter estimation. */
class FragmentModelInterval
{
    public:
        enum IntervalType
        {
            INTERGENIC,
            EXONIC,
            THREE_PRIME_EXONIC
        } type;

        FragmentModelInterval(const Interval& interval,
                              IntervalType type)
            : type(type)
            , seqname(interval.seqname)
            , start(interval.start)
            , end(interval.end)
            , strand(interval.strand)
            , three_prime_dist(-1)
            , tlen(-1)
            , tid(-1)
        {
        }

        FragmentModelInterval(SeqName seqname,
                              pos_t start,
                              pos_t end,
                              strand_t strand,
                              pos_t three_prime_dist,
                              pos_t tlen,
                              IntervalType type)
            : type(type)
            , seqname(seqname)
            , start(start)
            , end(end)
            , strand(strand)
            , three_prime_dist(three_prime_dist)
            , tlen(tlen)
            , tid(-1)
        {
        }
                               

        void add_alignment(const bam1_t* b)
        {
            rs.add_alignment(b);
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

        // Exon's distance from the 3' end of the transcript.
        pos_t three_prime_dist;

        // Transcript length
        pos_t tlen;

    private:
        /* A the sequence ID assigned by the BAM file, so we can arrange
         * intervals in the same order. */
        int32_t tid;

        friend void sam_scan(std::vector<FragmentModelInterval*>& intervals,
                             AlnCountTrie& T,
                             PosTable& mate1_pos_tab,
                             PosTable& mate2_pos_tab,
                             Queue<FragmentModelInterval*>& q,
                             const char* bam_fn,
                             const char* task_name);
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
              AlnCountTrie& T,
              PosTable& mate1_pos_tab,
              PosTable& mate2_pos_tab,
              Queue<FragmentModelInterval*>& q,
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

    sort(intervals.begin(), intervals.end(), FragmentModelIntervalPtrCmp());

    /* First interval which the current read may be contained in. */
    size_t j, j0 = 0;
    size_t n = intervals.size();

    size_t last_file_pos = 0, file_pos;
    size_t read_num = 0;

    /* Read the reads. */
    bam1_t* b = bam_init1();
    int32_t last_tid = -1;
    int32_t last_pos = -1;
    while (samread(bam_f, b) > 0) {
        ++read_num;

        if (read_num % 1000 == 0) {
            file_pos = samtell(bam_f);
            if (file_pos >= last_file_pos + input_block_size && input_size > 0) {
                Logger::get_task(task_name).inc();
                last_file_pos = file_pos;
            }
        }

        if (b->core.flag & BAM_FUNMAP || b->core.tid < 0) continue;

        if (b->core.tid < last_tid ||
            (b->core.tid == last_tid && b->core.pos < last_pos)) {
            Logger::abort(
                    "Excuse me, but I must insist that your SAM/BAM file be sorted. "
                    "Please run: 'samtools sort'.");
        }
        last_tid = b->core.tid;
        last_pos = b->core.pos;

        /* Count numbers of alignments by read. */
        if (b->core.flag & BAM_FREAD2) {
            if (b->core.n_cigar == 1) mate2_pos_tab.add(b, bam_f);
            T.inc_mate2(bam1_qname(b));
        }
        else {
            if (b->core.n_cigar == 1) mate1_pos_tab.add(b, bam_f);
            T.inc_mate1(bam1_qname(b));
        }

        /* Add reads to intervals in which they are contained. */
        for (j = j0; j < n; ++j) {
            if (b->core.tid < intervals[j]->tid) break;
            if (b->core.tid > intervals[j]->tid) {
                assert(j == j0);
                q.push(intervals[j0++]);
                continue;
            }

            if (b->core.pos < intervals[j]->start) break;
            if (b->core.pos > intervals[j]->end) {
                if (j == j0) {
                    q.push(intervals[j0++]);
                }
                continue;
            }

            pos_t b_end = (pos_t) bam_calend(&b->core, bam1_cigar(b)) - 1;
            if (b_end <= intervals[j]->end) {
                intervals[j]->add_alignment(b);
            }
        }
    }

    for(; j0 < n; ++j0) q.push(intervals[j0]);

    bam_destroy1(b);
    samclose(bam_f);

    if (task_name) Logger::pop_task(task_name);
}


/* A thread used to accumulate statistics to estimate the parameters of
 * the fragment model. */
class FragmentModelThread
{
    public:
        FragmentModelThread(Queue<FragmentModelInterval*>& q)
            : three_prime_dist_counts_0(
                    constants::transcript_3p_dist_len +
                    constants::transcript_3p_dist_pad, 0)
            , three_prime_dist_counts_1(
                    constants::transcript_3p_dist_len +
                    constants::transcript_3p_dist_pad, 0)
            , q(q)
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
                else if (interval->type == FragmentModelInterval::EXONIC) {
                    interval->rs.make_unique_read_counts(read_counts);
                    measure_fragment_lengths(read_counts);
                    measure_strand_bias(interval->strand, read_counts);
                    read_counts.clear();
                }
                else if (interval->type == FragmentModelInterval::THREE_PRIME_EXONIC) {
                    interval->rs.make_unique_read_counts(read_counts);
                    measure_three_prime_dist(interval->start,
                                             interval->end,
                                             interval->strand,
                                             interval->three_prime_dist,
                                             interval->tlen,
                                             read_counts);
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

        /* Counts of fragment starts at the given distance from the annotated
         * tss/tts. */
        std::vector<unsigned int> three_prime_dist_counts_0;
        std::vector<unsigned int> three_prime_dist_counts_1;

    private:
        void measure_fragment_lengths(const ReadSet::UniqueReadCounts& counts)
        {
            boost::unordered_map<unsigned int, unsigned int>::iterator c;
            ReadSet::UniqueReadCounts::const_iterator i;
            for (i = counts.begin(); i != counts.end(); ++i) {
                AlignedReadIterator j(*i->first);
                for (; j != AlignedReadIterator(); ++j) {
                    if (j->mate1 == NULL || j->mate2 == NULL) continue;
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

        void measure_three_prime_dist(
                pos_t start, pos_t end, strand_t strand, pos_t three_prime_dist,
                pos_t tlen, const ReadSet::UniqueReadCounts& counts)
        {
            ReadSet::UniqueReadCounts::const_iterator i;
            for (i = counts.begin(); i != counts.end(); ++i) {
                AlignedReadIterator j(*i->first);
                for (; j != AlignedReadIterator(); ++j) {
                    if (!j->mate1) continue;
                    pos_t d = j->mate1->strand == strand_pos ?
                                  j->mate1->start : j->mate1->end;
                    if (strand == strand_pos) {
                        d = three_prime_dist + end - d;
                    }
                    else {
                        d = three_prime_dist + d - start;
                    }

                    d = d * constants::transcript_3p_dist_len / tlen;

                    if (0 <= d && d < constants::transcript_3p_dist_len +
                                      constants::transcript_3p_dist_pad) {
                        if (strand == j->mate1->strand) {
                            ++three_prime_dist_counts_0[d];
                        }
                        else {
                            ++three_prime_dist_counts_1[d];
                        }
                    }
                }
            }
        }

        Queue<FragmentModelInterval*>& q;
        boost::thread* t;
};


FragmentModel::FragmentModel()
    : sb(NULL)
    , frag_len_dist(NULL)
{
    three_prime_dist[0] = three_prime_dist[1] = NULL;
}


FragmentModel::~FragmentModel()
{
    if (sb) delete sb;
    if (frag_len_dist) delete frag_len_dist;
    if (three_prime_dist[0]) delete three_prime_dist[0];
    if (three_prime_dist[1]) delete three_prime_dist[1];
}


void FragmentModel::estimate(TranscriptSet& ts,
                             const char* bam_fn,
                             const char* fa_fn)
{
    Queue<FragmentModelInterval*> q(constants::max_estimate_queue_size);

    std::vector<FragmentModelInterval*> intervals;

    std::vector<Interval> exonic;
    ts.get_consensus_exonic(exonic);
    std::vector<Interval>::iterator interval;
    for (interval = exonic.begin(); interval != exonic.end(); ++interval) {
        intervals.push_back(new FragmentModelInterval(
                                    *interval,
                                    FragmentModelInterval::EXONIC));
    }


    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        pos_t tlen = t->exonic_length();
        if (tlen >= constants::transcript_3p_dist_len +
                    constants::transcript_3p_dist_pad) {
            pos_t d = 0;
            if (t->strand == strand_neg) {
                for (Transcript::iterator e = t->begin(); e != t->end(); ++e) {
                    intervals.push_back(new FragmentModelInterval(
                                t->seqname,
                                e->start,
                                e->end,
                                t->strand,
                                d,
                                tlen,
                                FragmentModelInterval::THREE_PRIME_EXONIC));
                    d += e->end - e->start + 1;
                    if (d >= constants::transcript_3p_dist_len +
                             constants::transcript_3p_dist_pad) break;
                }
            }
            else {
                for (Transcript::reverse_iterator e = t->rbegin(); e != t->rend(); ++e) {
                    intervals.push_back(new FragmentModelInterval(
                                t->seqname,
                                e->start,
                                e->end,
                                t->strand,
                                d,
                                tlen,
                                FragmentModelInterval::THREE_PRIME_EXONIC));
                    d += e->end - e->start + 1;
                    if (d >= constants::transcript_3p_dist_len +
                             constants::transcript_3p_dist_pad) break;
                }
            }
        }
    }


    std::vector<FragmentModelThread*> threads;
    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads.push_back(new FragmentModelThread(q));
        threads.back()->start();
    }

    AlnCountTrie* alncnt = new AlnCountTrie();
    PosTable mate1_pos_tab, mate2_pos_tab;
    sam_scan(intervals, *alncnt,
             mate1_pos_tab, mate2_pos_tab, q,
             bam_fn, "Indexing reads");

    /* Index multireads. */
    unsigned long total_reads = alncnt->size();
    for (AlnCountTrieIterator i(*alncnt); i != AlnCountTrieIterator(); ++i) {
        if (i->second.first > constants::max_alignments ||
            i->second.second > constants::max_alignments) {
            blacklist.add(i->first);
        }
        else if (i->second.first > 1 || i->second.second > 1) {
            multireads.add(i->first);
        }
    }
    delete alncnt;

    Logger::info("Reads: %lu, %0.1f%% with multiple alignments",
                 total_reads,
                 100.0 * (double) multireads.size() / (double) total_reads);

    if (fa_fn != NULL) {
        sb = new sequencing_bias(fa_fn,
                                 mate1_pos_tab, mate2_pos_tab,
                                 constants::seqbias_num_reads,
                                 constants::seqbias_left_pos,
                                 constants::seqbias_right_pos,
                                 1.0);
    }
    else sb = NULL;

    for (size_t i = 0; i < constants::num_threads; ++i) {
        q.push(NULL);
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        threads[i]->join();
    }

    /* Collect 3' dist distribution statistics */
    std::vector<unsigned int> three_prime_dist_vals_0(
            constants::transcript_3p_dist_len + constants::transcript_3p_dist_pad);
    std::vector<unsigned int> three_prime_dist_counts_0(
            constants::transcript_3p_dist_len + constants::transcript_3p_dist_pad, 0);

    std::vector<unsigned int> three_prime_dist_vals_1(
            constants::transcript_3p_dist_len + constants::transcript_3p_dist_pad);
    std::vector<unsigned int> three_prime_dist_counts_1(
            constants::transcript_3p_dist_len + constants::transcript_3p_dist_pad, 0);

    for (unsigned int i = 0; i < three_prime_dist_vals_0.size(); ++i) {
        three_prime_dist_vals_0[i] = i;
        three_prime_dist_vals_1[i] = i;
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        for (size_t j = 0; j < threads[i]->three_prime_dist_counts_0.size(); ++j) {
            three_prime_dist_counts_0[j] += threads[i]->three_prime_dist_counts_0[j];
            three_prime_dist_counts_1[j] += threads[i]->three_prime_dist_counts_1[j];
        }
    }

    three_prime_dist[0] = new EmpDist(&three_prime_dist_vals_0.at(0),
                                      &three_prime_dist_counts_0.at(0),
                                      three_prime_dist_vals_0.size(),
                                      constants::transcript_3p_dist_w);
    three_prime_dist[1] = new EmpDist(&three_prime_dist_vals_1.at(0),
                                      &three_prime_dist_counts_1.at(0),
                                      three_prime_dist_vals_0.size(),
                                      constants::transcript_3p_dist_w);


    // Fit a line to the 3' bias to extrapolate past the end.
    size_t n = constants::transcript_3p_dist_len - constants::transcript_3p_dist_pad;
    std::vector<double> three_prime_dist_vals_d(n);
    std::vector<double> three_prime_dist_weight_d(n);
    for (unsigned int i = 0; i < n; ++i) {
        three_prime_dist_vals_d[i] = (double) (i + constants::transcript_3p_dist_pad);
    }

    double cov00, cov01, cov10, sumsq;
    for (unsigned int i = 0; i < n; ++i) {
        three_prime_dist_weight_d[i] = three_prime_dist[0]->pdf(
                i + constants::transcript_3p_dist_pad);
    }
    gsl_fit_linear(&three_prime_dist_vals_d.at(0), 1,
                   &three_prime_dist_weight_d.at(0), 1,
                   n, &three_prime_dist_c0[0], &three_prime_dist_c1[0],
                   &cov00, &cov01, &cov10, &sumsq);
    Logger::info("3' + : c0 = %e, c1 = %e",
                 three_prime_dist_c0[0],
                 three_prime_dist_c1[0]);


    for (unsigned int i = 0; i < three_prime_dist_vals_d.size(); ++i) {
        three_prime_dist_weight_d[i] = three_prime_dist[1]->pdf(
                i + constants::transcript_3p_dist_pad);
    }
    gsl_fit_linear(&three_prime_dist_vals_d.at(0), 1,
                   &three_prime_dist_weight_d.at(0), 1,
                   n, &three_prime_dist_c0[1], &three_prime_dist_c1[1],
                   &cov00, &cov01, &cov10, &sumsq);
    Logger::info("3' - : c0 = %e, c1 = %e",
                 three_prime_dist_c0[1],
                 three_prime_dist_c1[1]);


    // XXX: Debugging
    {
        FILE* out = fopen("3p_dist.0.tsv", "w");
        for (int i = 0; i < constants::transcript_3p_dist_len; ++i) {
            fprintf(out, "%d\t%f\n", i, three_prime_dist[0]->pdf(i));
        }
        fclose(out);

        out = fopen("3p_dist.1.tsv", "w");
        for (int i = 0; i < constants::transcript_3p_dist_len; ++i) {
            fprintf(out, "%d\t%f\n", i, three_prime_dist[1]->pdf(i));
        }
        fclose(out);
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
    Logger::info("Strand specificity: %0.1f%%", 100.0 * (1.0 + negentropy));

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

        /* prune fragment lengths with very low probability */
        for (i = 0; i < frag_lens.size(); ++i) {
            if ((float) frag_len_lens[i] / (float) sum_lens <
                    constants::min_frag_len_pr) {
                frag_len_lens[i] = 0;
            }
        }

        frag_len_dist = new EmpDist(frag_len_vals, frag_len_lens, frag_lens.size(),
                                    constants::frag_len_dist_smoothing);

        delete [] frag_len_vals;
        delete [] frag_len_lens;

        Logger::info("Fragment length distribution estimated.");
        Logger::info("Median fragment-length: %d", (int) frag_len_dist->median());
    }
    else {
        Logger::info("Too few paired-end reads to estimate fragment length distribution.");
    }

    for (size_t i = 0; i < constants::num_threads; ++i) {
        delete threads[i];
    }
}


float FragmentModel::frag_len_p(pos_t frag_len)
{
    if (frag_len_dist) return frag_len_dist->pdf(frag_len);
    else return gsl_ran_gaussian_pdf((double) frag_len - constants::frag_len_mu,
                                     constants::frag_len_sd);
}


float FragmentModel::frag_len_c(pos_t frag_len)
{
    if (frag_len_dist) return frag_len_dist->cdf(frag_len);
    else return gsl_cdf_gaussian_P((double) frag_len - constants::frag_len_mu,
                                   constants::frag_len_sd);
}


float FragmentModel::frag_len_med()
{
    if (frag_len_dist) return frag_len_dist->median();
    else return constants::frag_len_mu;
}


