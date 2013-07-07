
#include "pos_table.h"



#define NUM_PRIMES 28
static const uint32_t primes[NUM_PRIMES] = {
           53U,         97U,        193U,        389U,
          769U,       1543U,       3079U,       6151U,
        12289U,      24593U,      49157U,      98317U,
       196613U,     393241U,     786433U,    1572869U,
      3145739U,    6291469U,   12582917U,   25165843U,
     50331653U,  100663319U,  201326611U,  402653189U,
    805306457U, 1610612741U, 3221225473U, 4294967291U };


static const double max_load = 0.75;

/* marks a vacant cell */
static const int32_t nilpos = -1;


#ifndef MIN
#define MIN(a,b) ((a)<=(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>=(b)?(a):(b))
#endif


/* From Thomas Wang (http://www.cris.com/~Ttwang/tech/inthash.htm) */
static uint32_t hash(uint32_t a)
{
    a = (a ^ 61) ^ (a >> 16);
    a = a + (a << 3);
    a = a ^ (a >> 4);
    a = a * 0x27d4eb2d;
    a = a ^ (a >> 15);
    return a;
}



/* simple quadratic probing */
static uint32_t probe(uint32_t h, uint32_t i)
{
    const double c1 = 0.5;
    const double c2 = 0.5;

    return h
           + (uint32_t)(c1 * (double)i)
           + (uint32_t)(c2 * (double)(i*i));
}


static void pos_subtable_create(struct pos_subtable* T)
{
    T->n = 0;
    T->m = 0;
    T->A = malloc(sizeof(struct pos_table_val)*primes[T->n]);
    size_t i;
    for (i = 0; i < primes[T->n]; i++) {
        T->A[i].pos = nilpos;
        T->A[i].count = 0;
    }
    T->max_m = (size_t)(((double)primes[T->n]) * max_load);
}


static void pos_subtable_copy(struct pos_subtable* T, const struct pos_subtable* U)
{
    T->n     = U->n;
    T->m     = U->m;
    T->max_m = U->max_m;
    T->A     = malloc(sizeof(struct pos_table_val)*primes[T->n]);

    size_t i;
    for (i = 0; i < primes[T->n]; i++) {
        T->A[i].pos   = U->A[i].pos;
        T->A[i].count = U->A[i].pos;
    }
}


static void pos_subtable_destroy(struct pos_subtable* T)
{
    free(T->A);
    T->A = NULL;
}


static void pos_subtable_rehash(struct pos_subtable* T, size_t new_n);


static bool pos_subtable_inc(struct pos_subtable* T, int32_t pos)
{
    if (T->m == T->max_m) pos_subtable_rehash(T, T->n + 1);

    uint32_t h = hash((uint32_t)pos);
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while (T->A[j].pos != nilpos && T->A[j].pos != pos) {
        j = probe(h, ++i) % primes[T->n];
    }

    if (T->A[j].pos == nilpos) {
        T->A[j].pos = pos;
        T->A[j].count = 1;
        T->m++;
        return true;
    }
    else {
        T->A[j].count++;
        return false;
    }
}


static void pos_subtable_set(struct pos_subtable* T, int32_t pos, uint32_t count)
{
    if (T->m == T->max_m) pos_subtable_rehash(T, T->n + 1);

    uint32_t h = hash((uint32_t)pos);
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while (T->A[j].pos != nilpos && T->A[j].pos != pos) {
        j = probe(h, ++i) % primes[T->n];
    }

    if (T->A[j].pos == nilpos) {
        T->A[j].pos = pos;
        T->A[j].count = count;
    }
    else {
        T->A[j].count = count;
    }
}


static void pos_subtable_rehash(struct pos_subtable* T, size_t new_n)
{
    if (new_n >= NUM_PRIMES) {
        fprintf(stderr, "A table has grown too large!\n");
        exit(EXIT_FAILURE);
    }

    struct pos_subtable U;
    U.n = new_n;
    U.A = malloc(sizeof(struct pos_table_val) * primes[U.n]);
    size_t i;
    for (i = 0; i < primes[U.n]; i++) {
        U.A[i].pos = nilpos;
        U.A[i].count = 0;
    }

    U.m = 0;
    U.max_m = (size_t)(((double)primes[U.n]) * max_load);


    for (i = 0; i < primes[T->n]; i++) {
        if (T->A[i].pos == nilpos) continue;
        pos_subtable_set(&U, T->A[i].pos, T->A[i].count);
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = U.max_m;
}



static uint32_t pos_subtable_count(struct pos_subtable* T, int32_t pos)
{
    uint32_t h = hash((uint32_t)pos);
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while (T->A[j].pos != nilpos && T->A[j].pos != pos) {
        j = probe(h, ++i) % primes[T->n];
    }

    if (T->A[j].pos == pos) return T->A[j].count;
    else                     return 0;
}


void pos_table_create(struct pos_table* T, size_t n)
{
    T->seq_names = NULL;
    T->n = n;
    T->m = 0;

    T->ts[0] = malloc(n * sizeof(struct pos_subtable));
    T->ts[1] = malloc(n * sizeof(struct pos_subtable));

    size_t i, j;
    for (i = 0; i <= 1; i++) {
        for (j = 0; j < n; j++) {
            pos_subtable_create(&T->ts[i][j]);
        }
    }
}

void pos_table_copy(struct pos_table* T, const struct pos_table* U)
{
    T->seq_names = U->seq_names;
    T->n         = U->n;
    T->m         = U->m;

    T->ts[0] = malloc(T->n * sizeof(struct pos_subtable));
    T->ts[1] = malloc(T->n * sizeof(struct pos_subtable));

    size_t i, j;
    for (i = 0; i <= 1; i++) {
        for (j = 0; j < T->n; j++) {
            pos_subtable_copy(&T->ts[i][j], &U->ts[i][j]);
        }
    }
}


void pos_table_destroy(struct pos_table* T)
{
    size_t i, j;
    for (i = 0; i <= 1; i++) {
        for (j = 0; j < T->n; j++) {
            pos_subtable_destroy(&T->ts[i][j]);
        }
    }

    free(T->ts[0]);
    free(T->ts[1]);

    T->n = 0;
}


void pos_table_inc(struct pos_table* T, bam1_t* read)
{
    int32_t pos;
    if (bam1_strand(read)) pos = bam_calend2(&read->core, bam1_cigar(read)) - 1;
    else                   pos = read->core.pos;

    pos_table_inc_pos(T, read->core.tid, pos, bam1_strand(read));
}




void pos_table_inc_pos(struct pos_table* T, int32_t tid, int32_t pos, uint32_t strand)
{
    if (tid < 0 || (size_t) tid >= T->n) return;
    if (pos_subtable_inc(&T->ts[strand][tid], pos)) T->m++;
}


uint32_t table_count(struct pos_table* T, bam1_t* read)
{
    int32_t pos;
    if (bam1_strand(read)) pos = bam_calend2(&read->core, bam1_cigar(read)) - 1;
    else                    pos = read->core.pos;

    return pos_table_count_pos(T, read->core.tid, pos, bam1_strand(read));
}


uint32_t pos_table_count_pos(struct pos_table* T, int32_t tid, int32_t pos, uint32_t strand)
{
    if (tid < 0 || (size_t) tid >= T->n) return 0;
    return pos_subtable_count(&T->ts[strand][tid], pos);
}



static int pos_table_val_compare(const void* p1, const void* p2)
{
    return ((struct pos_table_val*)p1)->pos - ((struct pos_table_val*)p2)->pos;
}


void read_counts_create(struct read_counts* C, const struct pos_table* T)
{
    C->n = T->n;
    C->m = T->m;
    C->seq_names   = T->seq_names;

    C->mss[0] = malloc(C->n * sizeof(size_t));
    C->mss[1] = malloc(C->n * sizeof(size_t));

    C->xss[0] = malloc(C->n * sizeof(struct pos_table_val*));
    C->xss[1] = malloc(C->n * sizeof(struct pos_table_val*));

    int32_t  tid;
    uint32_t strand;
    size_t i, j;

    size_t m, n;
    struct pos_table_val* ys;
    struct pos_table_val* xs;


    for (strand = 0; strand <= 1; strand++) {
        for (tid = 0; (size_t) tid < T->n; tid++) {
            m  = T->ts[strand][tid].m;
            n  = T->ts[strand][tid].n;
            ys = T->ts[strand][tid].A;
            xs = malloc(m * sizeof(struct pos_table_val));

            for (i = 0, j = 0; j < primes[n]; j++) {
                if (ys[j].pos != nilpos) {
                    xs[i].pos   = ys[j].pos;
                    xs[i].count = ys[j].count;
                    i++;
                }
            }

            qsort(xs, m, sizeof(struct pos_table_val), pos_table_val_compare);

            C->mss[strand][tid] = m;
            C->xss[strand][tid] = xs;
        }
    }
}

void read_counts_copy(struct read_counts* C, const struct read_counts* B)
{
    C->n = B->n;
    C->m = B->m;
    C->seq_names = B->seq_names;

    int32_t  tid;
    uint32_t strand;
    size_t siz;

    for (strand = 0; strand <= 1; strand++) {
        C->mss[strand] = malloc(C->n * sizeof(size_t));
        C->xss[strand] = malloc(C->n * sizeof(struct pos_table_val*));
        for (tid = 0; (size_t) tid < C->n; tid++) {
            C->mss[strand][tid] = B->mss[strand][tid];
            siz = C->mss[strand][tid] * sizeof(struct pos_table_val);
            C->xss[strand][tid] = malloc(siz);
            memcpy(C->xss[strand][tid], B->xss[strand][tid], siz);
        }
    }
}


void read_counts_destroy(struct read_counts* C)
{
    int32_t  tid;
    uint32_t strand;

    for (strand = 0; strand <= 1; strand++) {
        for (tid = 0; (size_t) tid < C->n; tid++) {
            free(C->xss[strand][tid]);
            C->xss[strand][tid] = NULL;
        }
    }

    free(C->mss[0]); C->mss[0] = NULL;
    free(C->mss[1]); C->mss[1] = NULL;

    free(C->xss[0]); C->xss[0] = NULL;
    free(C->xss[1]); C->xss[1] = NULL;
}



/* find an index i, such that
 *      xs[i-1].pos < start <= xs[i].pos
 * using binary search.
 */
size_t bisect(struct pos_table_val* xs, size_t m, int32_t start)
{
    size_t a = 0;
    size_t b = m;
    size_t i = 0;

    while (a <= b) {
        i = a + (b - a) / 2;

        if (xs[i].pos < start)                  a = i + 1;
        else if (i > 0 && start <= xs[i-1].pos) b = i - 1;

        /* xs[i-1].pos <= start <= xs[i] */
        else break;
    }

    return i;
}


void read_counts_count(const struct read_counts* C,
                        int32_t tid, int32_t start, int32_t end, uint32_t strand,
                        unsigned int* ys)
{
    struct pos_table_val* xs = C->xss[strand][tid];
    size_t               m  = C->mss[strand][tid];

    if (m == 0) return;

    size_t i = bisect(xs, m, start);

    memset(ys, 0, m*sizeof(unsigned int));
    while (i < m && xs[i].pos <= end) {
        ys[ xs[i].pos - start ] = xs[i].count;
        i++;
    }
}

unsigned int read_counts_total(const struct read_counts* C,
                                int32_t tid, int32_t start, int32_t end, uint32_t strand)
{
    struct pos_table_val* xs = C->xss[strand][tid];
    size_t               m  = C->mss[strand][tid];

    if (m == 0) return 0;

    size_t i = bisect(xs, m, start);

    unsigned int total = 0;
    while (i < m && xs[i].pos <= end) {
        total += xs[i].count;
        i++;
    }

    return total;
}


void read_count_occurances(const struct read_counts* C,
                            int32_t tid,  int32_t start, int32_t end, uint32_t strand,
                            uint64_t* ks, size_t max_k)
{
    struct pos_table_val* xs = C->xss[strand][tid];
    size_t                m  = C->mss[strand][tid];

    if (m == 0) return;

    size_t i;

    i = bisect(xs, m, start);
    uint64_t nonzeros = 0;

    while (i < m && xs[i].pos <= end) {
        if (xs[i].count <= max_k) ks[xs[i].count]++;
        nonzeros++;
        i++;
    }

    uint64_t zeros = (end - start + 1) - nonzeros;

    /* Ignore and leading or trailing zeros if we are at the start or end of the
     * sequence. Many genome assemblies have several kilobases of N's at the
     * beginning and end. Considering these will lead to deflated statistics. */
    if (start <= xs[0].pos) {
        zeros -= MIN(end, xs[0].pos) - start + 1;
    }

    if (end >= xs[m-1].pos) {
        zeros -= end - MAX(start, xs[m-1].pos) + 1;
    }

    ks[0] += zeros;
}



void pos_table_dump(struct pos_table* T, struct read_pos** A_, size_t* N_, size_t limit)
{
    struct read_pos* A;
    size_t N = 0;
    size_t i, j;


    for (i = 0; i <= 1; i++) {
        for (j = 0; j < T->n; j++) {
            N += T->ts[i][j].m;
        }
    }

    if (limit > 0 && N > limit) N = limit;

    A = malloc(N * sizeof(struct read_pos));


    size_t u = 0;
    size_t v;
    for (i = 0; i <= 1; i++) {
        for (j = 0; j < T->n; j++) {
            for (v = 0; v < primes[T->ts[i][j].n]; v++) {
                if (T->ts[i][j].A[v].pos != -1) {
                    A[u].tid = j;
                    A[u].strand = i;
                    A[u].pos    = T->ts[i][j].A[v].pos;
                    A[u].count  = T->ts[i][j].A[v].count;
                    u++;
                    if (u >= N) goto table_dump_finish;
                }
            }
        }
    }
table_dump_finish:

    *A_ = A;
    *N_ = u;
}



