
#ifndef ISOLATOR_COMMON_HPP
#define ISOLATOR_COMMON_HPP

/* A genomic position. */
typedef long pos_t;

/* Strandedness, or lack therof. */
/* TODO: stop being a bitch and make these capital. */
typedef enum {
    strand_pos = 0,
    strand_neg = 1,
    strand_na  = 2
} strand_t;

/* Like strcmp but impart a more natural ordering on sequence names. */
int seqname_compare(const char* u, const char* v);


#endif

