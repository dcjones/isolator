
#ifndef ISOLATOR_SAMTOOLS_EXTRA_H
#define ISOLATOR_SAMTOOLS_EXTRA_H

#ifdef __cplusplus
extern "C" {
#endif

/* This file is to make a couple tweaks to samtools without actually modifying
 * the source code. I need to expose one function and add a new one.
 */

#include <stdlib.h>

#include "faidx.h"
#include "sam.h"

/* Report offset into a sam/bam file. */
size_t samtell(samfile_t* fp);


/* exposing a function defined in samtools/bam_aux.c */
void bam_init_header_hash(bam_header_t *header);


/* Change the behavior of the faidx_fetch_seq function to be more useful. If
 * coordinates are outside the actual sequence, write N's, rather than adjusting
 * the start,end. */
char* faidx_fetch_seq_forced_lower(const faidx_t* fai, const char *c_name, int p_beg_i, int p_end_i);


/* Calculate the true read-start. The start position stored in SAM/BAM files
 * doesn't include leading non-match cigar ops, which is fucking stupid, so we
 * work around it here. */
uint32_t bam_truepos(const bam1_core_t *c, const uint32_t* cigar);


/* Calculate read end-point include soft-clipping. */
uint32_t bam_trueend(const bam1_core_t *c, const uint32_t* cigar);

#ifdef __cplusplus
}
#endif

#endif


