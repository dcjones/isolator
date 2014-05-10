
#include "samtools_extra.h"
#include "khash.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 16384)

#ifndef _NO_RAZF
#include "razf.h"
#else
#ifdef _WIN32
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif
#define RAZF FILE
#define razf_read(fp, buf, size) fread(buf, 1, size, fp)
#define razf_open(fn, mode) fopen(fn, mode)
#define razf_close(fp) fclose(fp)
#define razf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define razf_tell(fp) ftello(fp)
#endif
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

#include <ctype.h>


typedef struct {
	int32_t line_len, line_blen;
	int64_t len;
	uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)


struct __faidx_t {
	RAZF *rz;
	int n, m;
	char **name;
	khash_t(s) *hash;
};


char* faidx_fetch_seq_forced_lower(const faidx_t* fai, const char *c_name,
                                   int p_beg_i, int p_end_i )
{
	int l;
	char c;
    khiter_t iter;
    faidx1_t val;
    char* seq0;
    char* seq = NULL;

    iter = kh_get(s, fai->hash, c_name);
    if (iter == kh_end(fai->hash)) return 0;

    val = kh_value(fai->hash, iter);

    seq0 = seq = malloc((p_end_i - p_beg_i + 2) * sizeof(char));
    if (seq0 == NULL) {
        fprintf(stderr, "Out of memory. Could not allocate %d bytes.\n",
                p_end_i - p_beg_i + 2);
        exit(EXIT_FAILURE);
    }
    seq0[p_end_i-p_beg_i+1] = '\0';

    /* entirely off the map: all Ns */
    if (p_beg_i >= (int)val.len || p_end_i < 0) {
        while (p_beg_i <= p_end_i) {
            *seq++ ='n';
            ++p_beg_i;
        }
        return seq0;
    }

    /* beginning is off the map */
    while (p_beg_i < 0 && p_beg_i <= p_end_i) {
        *seq++ = 'n';
        p_beg_i++;
    }

    /* end is off the map */
    while (p_end_i >= (int)val.len) {
        seq[p_end_i-p_beg_i] = 'n';
        p_end_i--;
    }

    /* retrieve the sequence */
	l = 0;
	razf_seek(fai->rz, val.offset + p_beg_i / val.line_blen * val.line_len + p_beg_i % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < p_end_i - p_beg_i + 1) {
		if (isgraph(c)) seq[l++] = tolower(c);
    }

    while (p_beg_i + l <= p_end_i) seq[l++] = 'n';

    return seq0;
}


struct __tamFile_t {
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	uint64_t n_lines;
	int is_first;
};


size_t samtell(samfile_t* fp)
{
    if (fp->type & 1) { /* TYPE_BAM == 1 */
#ifdef _USE_KNETFILE
        return knet_tell(fp->x.bam->x.fpr);
#else
        ftello(fp->x.bam->file);
#endif
    }
    else {
        return (size_t) gzoffset(fp->x.tamr->fp);
    }
}


uint32_t bam_truepos(const bam1_core_t *c, const uint32_t* cigar)
{
    uint32_t pos = c->pos;
    uint32_t k;
    for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
        int len = cigar[k] >> BAM_CIGAR_SHIFT;
        if (op != BAM_CMATCH) {
            pos -= len;
        }
        else break;
    }
    return pos;
}


uint32_t bam_trueend(const bam1_core_t *c, const uint32_t* cigar)
{
	uint32_t k, end;
    end = bam_truepos(c, cigar);
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CSOFT_CLIP)
			end += cigar[k] >> BAM_CIGAR_SHIFT;
	}
	return end;
}


