/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#include "sequencing_bias.hpp"
#include "logger.hpp"
#include "twobitseq.hpp"
#include "samtools/faidx.h"
#include "samtools/samtools_extra.h"
#include <climits>
#include <cmath>
#include <cctype>
#include <ctime>
#include <fstream>
#include <algorithm>

using namespace std;


static char complement(char c)
{
    switch( c ) {
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'n': return 'n';
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default:  return 'n';
    }
}

static void seqrc(char* seq, int n)
{
    char c;
    int i,j;
    i = 0;
    j = n-1;
    while( i < j ) {
        c = complement(seq[i]);
        seq[i] = complement(seq[j]);
        seq[j] = c;
        i++; j--;
    }

    if( i == j ) seq[i] = complement(seq[i]);
}



#if 0
/* round away from zero */
static double round_away(double a)
{
    if (a < 0.0) return floor(a);
    else         return ceil(a);
}
#endif


#if 0
/* simple uniform random numbers */
static double rand_uniform(double a, double b)
{
    return a + b * (double)rand() / (double)RAND_MAX;
}
#endif

#if 0
/* random gaussians (copied from GSL, to avoid dependency) */
static double rand_gauss(double sigma)
{
    double x, y, r2;

    do
    {
        /* choose x,y in uniform square (-1,-1) to (+1,+1) */
        x = -1 + 2 * rand_uniform(0.0,1.0);
        y = -1 + 2 * rand_uniform(0.0,1.0);

        /* see if it is in the unit circle */
        r2 = x * x + y * y;
    }
    while (r2 > 1.0 || r2 == 0);

    /* Box-Muller transform */
    return sigma * y * sqrt (-2.0 * log (r2) / r2);
}
#endif


double gauss_pdf (const double x, const double sigma)
{
  double u = x / fabs (sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-u * u / 2);
  return p;
}


struct ReadPosSeqnameCmp
{
    bool operator () (const ReadPos& a, const ReadPos& b)
    {
        return a.seqname < b.seqname;
    }
};


struct ReadPosCountCmp
{
    bool operator () (const ReadPos& a, const ReadPos& b)
    {
        return a.count > b.count;
    }
};


/*
static int read_pos_tid_count_compare(const void* p1, const void* p2)
{
    int c = (int)(((read_pos*)p1)->tid - ((read_pos*)p2)->tid);

    if (c == 0) {
        return (int)((read_pos*)p2)->count - (int)((read_pos*)p1)->count;
    }
    else return c;
}
*/



sequencing_bias::sequencing_bias()
    : ref_f(NULL)
    , M(NULL)
{
}


sequencing_bias::sequencing_bias(const char* ref_fn,
                                 PosTable& T,
                                 size_t max_reads,
                                 pos_t L, pos_t R,
                                 const char* task_name,
                                 double complexity_penalty)
    : ref_f(NULL)
    , M(NULL)
{
    build(ref_fn, T, max_reads, L, R,
          task_name,
          complexity_penalty);
}


void sequencing_bias::clear()
{
    if (ref_f) {
        fai_destroy(ref_f);
        ref_f = NULL;
    }
    ref_fn.clear();

    delete M;
    M = NULL;
}



void sequencing_bias::build(const char* ref_fn,
                            PosTable& T,
                            size_t max_reads,
                            pos_t L, pos_t R,
                            const char* task_name,
                            double complexity_penalty)
{
    Logger::push_task(task_name);

    clear();
    const size_t min_positions = 1000;
    if (T.size() < min_positions) return;


    this->ref_fn = ref_fn;

    this->L = L;
    this->R = R;

    const size_t max_dump = 10000000;
    std::vector<ReadPos> S;
    S.reserve(max_dump);
    T.dump(S, max_dump);

    /* sort by tid (so we can load one chromosome at a time) */
    random_shuffle(S.begin(), S.end());
    sort(S.begin(), S.end(), ReadPosSeqnameCmp());

    //sort(S.begin(), S.end(), ReadPosCountCmp());
    //sort(S.begin(), S.begin() + max_reads, ReadPosSeqnameCmp());

    /* sample foreground and background kmer frequencies */
    ref_f = fai_load(ref_fn);
    if (ref_f == NULL) {
        Logger::abort("Can't open fasta file '%s'.", ref_fn);
    }

    std::deque<twobitseq*> foreground_seqs;
    std::deque<twobitseq*> background_seqs;

    /* background sampling */
    int bg_samples = 1; // make this many samples for each read
    int bg_sample_num; // keep track of the number of samples made
    pos_t bg_pos;

    int            seqlen    = 0;
    SeqName        curr_seqname;
    char*          seq       = NULL;

    char* local_seq;
    local_seq = new char[ L + R + 2 ];
    local_seq[L+R+1] = '\0';

    std::vector<ReadPos>::iterator i;
    for (i = S.begin(); i != S.end() && i != S.begin() + max_reads; ++i) {

        /* Load/switch sequences (chromosomes) as they are encountered in the
         * read stream. The idea here is to avoid thrashing by loading a large
         * sequence, but also avoid overloading memory by only loading one
         * chromosome at a time. */
        if (i->seqname != curr_seqname) {
            if (seq) free(seq);

            seq = faidx_fetch_seq(ref_f, i->seqname.get().c_str(), 0, INT_MAX, &seqlen);
            Logger::debug("read sequence %s.", i->seqname.get().c_str());

            if (seq == NULL) {
                Logger::warn("warning: reference sequence not found, skipping.");
            }
            else {
                for (char* c = seq; *c; c++) *c = tolower(*c);
            }

            curr_seqname = i->seqname;
        }

        if (seq == NULL) continue;

        /* add a foreground sequence */
        if (i->strand == strand_neg) {
            if (i->pos < R || i->pos >= seqlen - L) continue;
            memcpy(local_seq, seq + i->pos - R, (L+1+R)*sizeof(char));
            seqrc(local_seq, L+1+R);
        }
        else {
            if (i->pos < L || i->pos >= seqlen - R) continue;
            memcpy(local_seq, seq + (i->pos - L), (L+1+R)*sizeof(char));
        }

        if (strchr(local_seq, 'n') != NULL) continue;

        foreground_seqs.push_back(new twobitseq(local_seq));


        /* add a background sequence */
        /* adjust the current read position randomly, and sample */
        for (bg_sample_num = 0; bg_sample_num < bg_samples;) {
            bg_pos = random_uniform_int(rng,
                boost::random::uniform_int_distribution<pos_t>::param_type(i->start, i->end));

            if (i->strand == strand_neg) {
                if (bg_pos < R || bg_pos >= seqlen - L) continue;
                memcpy(local_seq, seq + bg_pos - R, (L+1+R)*sizeof(char));
                seqrc(local_seq, L+1+R);
            }
            else {
                if (bg_pos < L || bg_pos >= seqlen - R) continue;
                memcpy(local_seq, seq + (bg_pos-L), (L+1+R)*sizeof(char));
            }

            if (strchr(local_seq, 'n') != NULL) continue;

            background_seqs.push_back(new twobitseq(local_seq));
            bg_sample_num++;
        }
    }


    size_t max_parents  = 4;
    size_t max_distance = 10;

    /* A bit of a hack: if we are training on very few reads (a couple thousand,
     * as a opposed to tens of thousands), we tend to end up with too sparse of
     * a model. */
    if (foreground_seqs.size() < 10000) complexity_penalty = 0.25;

    M = new motif(background_seqs,
                  foreground_seqs,
                  L + 1 + R,
                  max_parents,
                  max_distance,
                  complexity_penalty,
                  task_name);

    std::deque<twobitseq*>::iterator seqit;
    for (seqit = background_seqs.begin(); seqit != background_seqs.end(); seqit++) {
        delete *seqit;
    }

    for (seqit = foreground_seqs.begin(); seqit != foreground_seqs.end(); seqit++) {
        delete *seqit;
    }

    free(seq);
    delete [] local_seq;

    Logger::pop_task(task_name);
}


sequencing_bias::~sequencing_bias()
{
    clear();
}


double sequencing_bias::get_bias(const twobitseq& seq, pos_t pos) const
{
    if (M == NULL || pos < L || (pos_t) seq.size() - pos <= R) return 1.0;

    return M->eval(seq, pos - L);
}


string sequencing_bias::model_graph() const
{
    return M->model_graph((int) L);
}


