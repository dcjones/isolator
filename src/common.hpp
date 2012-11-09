
#ifndef ISOLATOR_COMMON_HPP
#define ISOLATOR_COMMON_HPP

#include <boost/flyweight.hpp>
#include <cstdlib>
#include <string>

#include "logger.hpp"

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

/* Define types for string interning using boost::flyweight. */
struct gene_id_tag {};
typedef boost::flyweight<std::string, boost::flyweights::tag<gene_id_tag> > GeneID;

struct transcript_id_tag {};
typedef boost::flyweight<std::string, boost::flyweights::tag<transcript_id_tag> > TranscriptID;

struct seq_name_tag {};
typedef boost::flyweight<std::string, boost::flyweights::tag<transcript_id_tag> > SeqName;

#endif

