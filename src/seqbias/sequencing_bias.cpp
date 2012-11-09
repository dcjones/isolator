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



/* round away from zero */
static double round_away(double a)
{
    if (a < 0.0) return floor(a);
    else         return ceil(a);
}


/* simple uniform random numbers */
static double rand_uniform(double a, double b)
{
    return a + b * (double)rand() / (double)RAND_MAX;
}

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


double gauss_pdf (const double x, const double sigma)
{
  double u = x / fabs (sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-u * u / 2);
  return p;
}



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
    , M1(NULL)
    , M2(NULL)
{}


#if 0
sequencing_bias::sequencing_bias(const char* model_fn)
    : ref_f(NULL)
    , M1(NULL)
    , M2(NULL)
{
    std::ifstream fin;
    fin.open(model_fn);

    if (!fin) {
        logger::abort("Can't open file %s for reading.\n", model_fn);
    }

    YAML::Parser parser(fin);
    YAML::Node doc;
    parser.GetNextDocument(doc);

    doc["L"] >> L;
    doc["R"] >> R;

    const YAML::Node* node;

    if ((node = doc.FindValue("motif1")) != NULL) M1 = new motif(*node);
    if ((node = doc.FindValue("motif2")) != NULL) M2 = new motif(*node);
}
#endif

#if 0
sequencing_bias::sequencing_bias(const char* ref_fn,
                                 const char* model_fn)
    : ref_f(NULL)
    , M1(NULL)
    , M2(NULL)

{
    std::ifstream fin;
    fin.open(model_fn);

    YAML::Parser parser(fin);
    YAML::Node doc;
    parser.GetNextDocument(doc);

    doc["L"] >> L;
    doc["R"] >> R;

    const YAML::Node* node;

    if ((node = doc.FindValue("motif1")) != NULL) M1 = new motif(*node);
    if ((node = doc.FindValue("motif2")) != NULL) M2 = new motif(*node);


    if (ref_fn != NULL) {
        this->ref_fn = ref_fn;
        ref_f = fai_load(ref_fn);
        if (ref_f == NULL) {
            logger::abort("Can't open indexed FASTA file %s\n", ref_fn);
        }
    }
    else ref_f = NULL;
}
#endif


sequencing_bias::sequencing_bias(const char* ref_fn,
                                 PosTable& T1, PosTable& T2,
                                 size_t max_reads,
                                 pos_t L, pos_t R,
                                 double complexity_penalty)
    : ref_f(NULL)
    , M1(NULL)
    , M2(NULL)
{
    build(ref_fn, T1, T2, max_reads, L, R,
          complexity_penalty);
}




#if 0
void sequencing_bias::to_yaml(YAML::Emitter& out) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "L";
    out << YAML::Value << L;

    out << YAML::Key   << "R";
    out << YAML::Value << R;

    if (M1 != NULL) {
        out << YAML::Key   << "motif1";
        out << YAML::Value;
        M1->to_yaml(out);
    }

    if (M2 != NULL) {
        out << YAML::Key   << "motif2";
        out << YAML::Value;
        M2->to_yaml(out);
    }

    out << YAML::EndMap;
}
#endif


#if 0
void sequencing_bias::save_to_file(const char* fn) const
{
    FILE* f = fopen(fn, "w");
    if (f == NULL) {
        logger::abort("Can't open file %s for writing.", fn);
    }

    YAML::Emitter out;
    this->to_yaml(out);

    fputs(out.c_str(), f);
    fclose(f);
}
#endif


void sequencing_bias::clear()
{
    if (ref_f) {
        fai_destroy(ref_f);
        ref_f = NULL;
    }
    ref_fn.clear();

    delete M1;
    M1 = NULL;

    delete M2;
    M2 = NULL;
}



void sequencing_bias::build(const char* ref_fn,
                            PosTable& T1, PosTable& T2,
                            size_t max_reads,
                            pos_t L, pos_t R,
                            double complexity_penalty)
{
    clear();

    const size_t min_positions = 1000;

    if (T1.size() >= min_positions) {
        const char* task_name;
        if (T2.size() < min_positions) {
            task_name = "Estimating sequence bias";
        }
        else {
            task_name = "Estimating mate 1 sequence bias";
        }

        Logger::push_task(task_name);
        buildn(&M1, ref_fn, T1, max_reads, L, R, complexity_penalty, task_name);
        Logger::pop_task(task_name);
    }

    if (T2.size() >= min_positions) {
        const char* task_name = "Estimating mate 2 sequence bias";
        Logger::push_task(task_name);
        buildn(&M2, ref_fn, T2, max_reads, L, R, complexity_penalty, task_name);
        Logger::pop_task(task_name);
    }
}


struct ReadPosSeqnameCmp
{
    bool operator () (const ReadPos& a, const ReadPos& b)
    {
        return a.seqname < b.seqname;
    }
};


void sequencing_bias::buildn(motif** Mn,
                             const char* ref_fn,
                             PosTable& T,
                             size_t max_reads,
                             pos_t L, pos_t R,
                             double complexity_penalty,
                             const char* task_name)
{
    this->ref_fn = ref_fn;

    this->L = L;
    this->R = R;

    const size_t max_dump = 10000000;
    std::vector<ReadPos> S(max_dump);
    T.dump(S);

    /* sort by tid (so we can load one chromosome at a time) */
    random_shuffle(S.begin(), S.end());
    sort(S.begin(), S.end(), ReadPosSeqnameCmp());

    /* sample foreground and background kmer frequencies */
    ref_f = fai_load(ref_fn);
    if (ref_f == NULL) {
        Logger::abort("Can't open fasta file '%s'.", ref_fn);
    }

    std::deque<twobitseq*> foreground_seqs;
    std::deque<twobitseq*> background_seqs;


    /* background sampling */
    int bg_samples = 2; // make this many samples for each read
    int bg_sample_num;  // keep track of the number of samples made
    pos_t bg_pos;

    int            seqlen    = 0;
    SeqName        curr_seqname;
    char*          seq       = NULL;

    char* local_seq;
    local_seq = new char[ L + R + 2 ];
    local_seq[L+R+1] = '\0';

    std::vector<ReadPos>::iterator i;
    for (i = S.begin(); i != S.begin() + max_reads; ++i) {

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
            memcpy(local_seq, seq + (i->pos-L), (L+1+R)*sizeof(char));
        }

        if (strchr(local_seq, 'n') != NULL) continue;

        foreground_seqs.push_back(new twobitseq(local_seq));


        /* add a background sequence */
        /* adjust the current read position randomly, and sample */
        for (bg_sample_num = 0; bg_sample_num < bg_samples;) {

            bg_pos = i->pos + (pos_t)round_away(rand_gauss(500));

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

    *Mn = new motif(background_seqs,
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
    free(local_seq);
}


sequencing_bias::~sequencing_bias()
{
    clear();
}


double* sequencing_bias::get_bias(const char* seqname,
                                  pos_t start, pos_t end,
                                  strand_t strand) const
{
    return get_mate1_bias(seqname, start, end, strand);
}



double sequencing_bias::get_bias(const twobitseq& seq, pos_t pos) const
{
    return get_mate1_bias(seq, pos);
}



double* sequencing_bias::get_mate1_bias(const char* seqname,
                                        pos_t start, pos_t end,
                                        strand_t strand) const
{
    return get_maten_bias(M1, seqname, start, end, strand);
}


double sequencing_bias::get_mate1_bias(const twobitseq& seq, pos_t pos) const
{
    return get_maten_bias(M1, seq, pos);
}


double* sequencing_bias::get_mate2_bias(const char* seqname,
                                        pos_t start, pos_t end,
                                        strand_t strand) const
{
    return get_maten_bias(M2, seqname, start, end, strand);
}


double sequencing_bias::get_mate2_bias(const twobitseq& seq, pos_t pos) const
{
    return get_maten_bias(M2, seq, pos);
}


double* sequencing_bias::get_maten_bias(const motif* Mn,
                                        const char* seqname,
                                        pos_t start, pos_t end,
                                        strand_t strand) const
{
    if (strand == strand_na || ref_f == NULL || Mn == NULL) return NULL;

    pos_t i;
    pos_t seqlen = end - start + 1;

    double* bs = new double [seqlen];
    for (i = 0; i < seqlen; i++) bs[i] = 1.0;

    char* seqstr;

    if (strand == strand_neg) {
        seqstr = faidx_fetch_seq_forced_lower(ref_f, seqname,
                                              start - R, end + L);
        if (seqstr) seqrc(seqstr, seqlen + R + L);
    }
    else {
        seqstr = faidx_fetch_seq_forced_lower(ref_f, seqname,
                                              start - L, end + R);
    }


    if (!seqstr) return bs;

    twobitseq seq(seqstr);

    for (i = 0; i < seqlen; i++) {
        bs[i] = Mn->eval(seq, i);
    }


    free(seqstr);
    return bs;
}



double sequencing_bias::get_maten_bias(const motif* Mn,
                                       const twobitseq& seq, pos_t pos) const
{
    if (Mn == NULL || pos < L || (pos_t) seq.size() - pos <= R) return 1.0;
    return Mn->eval(seq, pos - L);
}




string sequencing_bias::model_graph() const
{
    return M1->model_graph((int) L);
    // TODO: also output M2
}


#if 0
kmer_matrix tabulate_bias(double* kl,
                          pos_t L, pos_t R, int k,
                          const char* ref_fn,
                          const char* reads_fn,
                          bool mate1,
                          bool mate2,
                          const char* model_fn)
{
    /* This function procedes very much like the sequencing_bias::build
     * function, but does not finally train a motif, but rather tabulates
     * k-mer frequencies. */
    size_t max_reads = 1000000;

    kmer_matrix dest((size_t) (L + 1 + R), (size_t) k);
    dest.set_all(0.0);

    faidx_t* ref_f = fai_load(ref_fn);
    if (ref_f == NULL) {
        Logger::abort("Can't open fasta file '%s'.", ref_fn);
    }

    samfile_t* reads_f = samopen(reads_fn, "rb", NULL);
    if (reads_f == NULL) {
        Logger::abort("Can't open bam file '%s'.", reads_fn);
    }

    sequencing_bias* sb = NULL;
    if (model_fn != NULL) {
        //sb = new sequencing_bias(ref_fn, model_fn);
        Logger::warn("Seqbias model serialization disabled.");
    }

    bam_init_header_hash(reads_f->header);
    bam1_t* read = bam_init1();

    pos_table T;
    pos_table_create(&T, reads_f->header->n_targets);
    T.seq_names = reads_f->header->target_name;

    size_t hashed_count = 0;
    while (samread(reads_f, read) > 0) {
        if (read->core.flag & BAM_FPAIRED) {
            if (read->core.flag & BAM_FREAD1 && !mate1) continue;
            if (read->core.flag & BAM_FREAD2 && !mate2) continue;
        }

        if (++hashed_count % 1000000 == 0) {
            Logger::info("hashed %zu reads.", hashed_count);
        }
        pos_table_inc(&T, read);
    }
    Logger::info("hashed %zu reads.", hashed_count);

    read_pos* S_tmp;
    read_pos* S;
    size_t N;
    const size_t max_dump = 10000000;
    pos_table_dump(&T, &S_tmp, &N, max_dump);

    /* sort by count */
    //qsort(S_tmp, N, sizeof(read_pos), read_pos_count_compare);

    /* consider only reads with at least one duplicate */
    size_t i;
    for (i = 0; i < N; ++i) {
        if (S_tmp[i].count <= 1) break;
    }

    /* (unless there are very few of these reads */
    if (i > 10000) {
        max_reads = std::min<size_t>(max_reads, i);
        Logger::info("%zu reads with duplicates.", i);
    }
    else {
        i = N;
    }

    /* ignore the top 1%, as they tend to be vastly higher than anything else,
     * and thus can throw things off when training with small numbers of reads
     * */
    S = S_tmp + i/100;
    max_reads = min<size_t>(max_reads, 99*i/100); 

    /* sort by tid (so we can load one chromosome at a time) */
    qsort(S, std::min<size_t>(max_reads, N), sizeof(read_pos), read_pos_tid_compare);

    twobitseq tbs;

    char* local_seq;
    local_seq = new char[ (k - 1) + L + 1 + R + 1 ];
    local_seq[(k - 1) + L + 1 + R] = '\0';

    double         w;
    char*          seq       = NULL;
    int            seqlen    = 0;
    int            curr_tid  = -1;
    char*          seqname   = NULL;

    for (i = 0; i < N && i < max_reads; i++) {
        if (S[i].tid != curr_tid) {
            seqname = T.seq_names[S[i].tid];
            free(seq);
            seq = faidx_fetch_seq(ref_f, seqname, 0, INT_MAX, &seqlen);
            Logger::info("read sequence %s.", seqname);
            curr_tid = S[i].tid;

            if (seq == NULL) {
                Logger::warn("warning: reference sequence not found, skipping.");
            }
        }

        if (seq == NULL) continue;

        if (S[i].strand == strand_neg) {
            if (S[i].pos < R || S[i].pos >= seqlen - L - (k-1)) continue;
            memcpy(local_seq, seq + S[i].pos - R, ((k-1)+L+1+R)*sizeof(char));
            seqrc(local_seq, L+1+R);
        }
        else {
            if (S[i].pos < L + (k-1) || S[i].pos >= seqlen - R) continue;
            memcpy(local_seq, seq + (S[i].pos - L - (k-1)), ((k-1)+L+1+R)*sizeof(char));
        }

        if (sb) {
            // TODO: get_bias
        }
        else w = 1.0;

        // XXX
        w = (double) S[i].count;

        tbs = local_seq;
        kmer K;
        for (pos_t pos = (k-1); pos < (k-1) + L + 1 + R; ++pos) {
            K = tbs.get_kmer(k, pos);
            dest(pos - (k-1), K) += w;
        }
    }

    /* compute KL divergence */

    /* estimate a background distribution by averaging across the entire window
     * */
    int four_to_k = 1 << (2 * k);
    double* bg = new double[four_to_k];
    memset(bg, 0, four_to_k * sizeof(double));
    for (pos_t pos = 0; pos < L + 1 + R; ++pos) {
        for (kmer K = 0; K < (kmer) four_to_k; ++K) {
            bg[K] += dest(pos, K);
        }
    }

    kmer_matrix norm_dest(dest);
    norm_dest.make_distribution();

    double z = 0.0;
    for (kmer K = 0; K < (kmer) four_to_k; ++K) z += bg[K];
    for (kmer K = 0; K < (kmer) four_to_k; ++K) bg[K] /= z;


    /* Compute the (symetric) kullback-leibler divegnnce */
    memset(kl, 0, (L + 1 + R) * sizeof(double));
    for (pos_t pos = 0; pos < L + 1 + R; ++pos) {
        kl[pos] = 0.0;
        for (kmer K = 0; K < (kmer) four_to_k; ++K) {
            if (norm_dest(pos, K) > 0.0) {
                kl[pos] += norm_dest(pos, K) * (log2(norm_dest(pos, K)) - log2(bg[K]));
            }

            if (bg[K] > 0.0) {
                kl[pos] += bg[K] * (log2(bg[K]) - log2(norm_dest(pos, K)));
            }
        }
    }

    delete [] bg;


    free(seq);
    free(local_seq);
    free(S_tmp);
    bam_destroy1(read);
    pos_table_destroy(&T);
    delete sb;
    samclose(reads_f);

    return dest;
}
#endif


