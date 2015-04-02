

/* A simple program to count reads. Similar to htseq-count, but much much
 * faster. */

#include <cstdio>
#include <getopt.h>

#include "alnindex.hpp"
#include "logger.hpp"
#include "read_set.hpp"
#include "transcripts.hpp"
#include "trie.hpp"
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

enum Stranded {
    STRANDED_NONE,
    STRANDED_FORWARD,
    STRANDED_REVERSE
};


static void print_usage(FILE* out)
{
    fprintf(out, "Usage: samcnt sam_file gtf_file\n");
}


static void print_help(FILE* out)
{
    print_usage(out);
    fprintf(out,
            "\n"
            "Arguments:\n"
            "  -h, --help                print this message\n"
            "  -o, --output=FILE         print the count table to this file\n"
            "  -s, --stranded=STRANDED   count in strand specific manner. STRANDED\n"
            "                            is one of \"yes\", \"no\", or \"reverse\".\n"
            "                            (default: no)\n"
            "  -t, --type=TYPE           GTF/GFF feature to use (default: exon)\n"
            "  -i, --id-attr=ID          sum counts that share this attribute\n"
            "                            (default: gene_id)\n"
            "  -I, --transcript-attr=ID  attribute used to define transcripts\n"
            "                            (default: transcript_id)\n"
            "\n");
}


struct MateCount {
    MateCount() : mate1_count(0), mate2_count(0) {}
    unsigned int mate1_count, mate2_count;
};


struct CountInterval {
    CountInterval(const TranscriptSetLocus& locus, int tid)
        : locus(locus)
        , tid(tid)
    {}

    bool operator < (const CountInterval& other) const
    {
        if (tid != other.tid) return tid < other.tid;
        else if (locus.min_start != other.locus.min_start) {
            return locus.min_start < other.locus.min_start;
        }
        else {
            return locus.max_end < other.locus.max_end;
        }
    }

    TranscriptSetLocus locus;
    int tid;
};


static void process_locus(TrieMap<unsigned long>& counts,
                          const TranscriptSetLocus& locus,
                          const ReadSet& rs, Stranded stranded)
{
    for (std::map<long, AlignedRead*>::const_iterator r = rs.rs.begin();
            r != rs.rs.end(); ++r) {
        AlignedReadIterator aln(*r->second);
        if (aln == AlignedReadIterator()) continue;
        const AlignmentPair& frag = *aln;
        GeneID gene_id;

        for (TranscriptSetLocus::const_iterator t = locus.begin();
             t != locus.end(); ++t) {

            if (stranded == STRANDED_FORWARD &&
                ((frag.mate1 && frag.mate1->strand != t->strand) ||
                 (frag.mate2 && frag.mate2->strand == t->strand))) {
                continue;
            }
            else if (stranded == STRANDED_REVERSE &&
                     ((frag.mate1 && frag.mate1->strand == t->strand) ||
                      (frag.mate2 && frag.mate2->strand != t->strand))) {
                continue;
            }

            if (frag.frag_len(*t) >= 0) {
                if (gene_id != GeneID() && gene_id != t->gene_id) {
                    gene_id = GeneID();
                    break;
                }
                else {
                    gene_id = t->gene_id;
                }
            }
        }

        if (gene_id != GeneID()) {
            counts[gene_id.get().c_str()] += 1;
        }
    }
}


int main(int argc, char* argv[])
{

    int opt, opt_idx;
    const char* output_filename = "-";

    static struct option long_options[] =
    {
        {"help",            no_argument,       NULL, 'h'},
        {"output",          required_argument, NULL, 'o'},
        {"stranded",        required_argument, NULL, 's'},
        {"type",            required_argument, NULL, 't'},
        {"transcript-attr", required_argument, NULL, 'I'},
        {"id-attr",         required_argument, NULL, 'i'},
        {0, 0, 0, 0}
    };

    const char* feature = "exon";
    const char* gid_attr = "gene_id";
    const char* tid_attr = "transcript_id";
    Stranded stranded = STRANDED_NONE;

    while (true) {
        opt = getopt_long(argc, argv, "ho:s:t:I:i:", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_help(stdout);
                return EXIT_FAILURE;

            case 'o':
                output_filename = optarg;
                break;

            case 's':
                if (strcmp(optarg, "yes") == 0 || strcmp(optarg, "forward") == 0) {
                    stranded = STRANDED_FORWARD;
                }
                else if (strcmp(optarg, "no") == 0 || strcmp(optarg, "none") == 0) {
                    stranded = STRANDED_NONE;
                }
                else if (strcmp(optarg, "reverse") == 0) {
                    stranded = STRANDED_REVERSE;
                }
                else {
                    fprintf(stderr, "\"%s\" is not a valid argument for \"--stranded\"\n", optarg);
                    return EXIT_FAILURE;
                }
                break;

            case 't':
                feature = optarg;
                break;

            case 'i':
                gid_attr = optarg;
                break;

            case 'I':
                tid_attr = optarg;
                break;

            case '?':
                fprintf(stderr, "\n");
                print_help(stderr);
                return EXIT_FAILURE;

            default:
                abort();
        }
    }

    FILE* output_file = stdout;
    if (strcmp(output_filename, "-") != 0) {
        output_file = fopen(output_filename, "w");
        if (!output_file) {
            Logger::abort("Unable to open %s for writing.", output_filename);
        }
    }

    if (optind + 1 >= argc) {
        print_usage(stderr);
        return EXIT_FAILURE;
    }

    const char* sam_filename = argv[optind];
    const char* gtf_filename = argv[optind+1];

    TranscriptSet transcripts;
    transcripts.read_gtf(gtf_filename, 0, false, feature, tid_attr, gid_attr);

    // Make one pass over the reads to count how many alignments each has.

    TrieMap<MateCount> aln_counts;

    samfile_t* sam_file;
    sam_file = samopen(sam_filename, "rb", NULL);
    if (!sam_file) {
        sam_file = samopen(sam_filename, "r", NULL);
    }
    if (!sam_file) {
        Logger::abort("Can't open SAM/BAM file %s.", sam_filename);
    }
    bam1_t* b = bam_init1();
    int last_tid = -1;
    pos_t last_pos = -1;
    while (samread(sam_file, b) >= 0) {
        if (b->core.tid == -1) continue;

        if (b->core.tid < last_tid ||
            (b->core.tid == last_tid && b->core.pos < last_pos)) {
            Logger::abort("Can't proceed. SAM/BAM file must be sorted.");
        }
        last_tid = b->core.tid;
        last_pos = b->core.pos;

        MateCount& count = aln_counts[bam1_qname(b)];
        if (b->core.flag & BAM_FREAD2) {
            ++count.mate2_count;
        }
        else {
            ++count.mate1_count;
        }
    }
    samclose(sam_file);

    // Make second pass to actually count reads.

    std::vector<CountInterval> intervals;
    for (TranscriptSetLocusIterator t(transcripts);
         t != TranscriptSetLocusIterator(); ++t) {
        intervals.push_back(CountInterval(*t, -1));
    }

    sam_file = samopen(sam_filename, "rb", NULL);
    if (!sam_file) {
        sam_file = samopen(sam_filename, "r", NULL);
    }
    if (!sam_file) {
        Logger::abort("Can't open SAM/BAM file %s.", sam_filename);
    }

    // sort intervals in the same order as the BAM file
    bam_init_header_hash(sam_file->header);
    khash_t(s)* tbl = reinterpret_cast<khash_t(s)*>(sam_file->header->hash);
    khiter_t k;

    for (std::vector<CountInterval>::iterator i = intervals.begin();
         i != intervals.end(); ++i) {
        const char* seqname = i->locus.seqname.get().c_str();
        k = kh_get(s, tbl, seqname);
        if (k == kh_end(tbl)) i->tid = -1;
        else                  i->tid = kh_value(tbl, k);
    }

    std::sort(intervals.begin(), intervals.end());

    TrieMap<unsigned long> counts;
    for (TranscriptSet::iterator t = transcripts.begin();
         t != transcripts.end(); ++t) {
        counts[t->gene_id.get().c_str()] = 0;
    }

    AlnIndex alnindex;
    ReadSet rs;
    size_t i = 0; // index into intervals

    while (samread(sam_file, b) >= 0 && i < intervals.size()) {
        pos_t start_pos = b->core.pos;
        pos_t end_pos = bam_calend(&b->core, bam1_cigar(b));

        long idx = alnindex.get(bam1_qname(b));

        // skip reads with multiple alignments
        MateCount& count = aln_counts[bam1_qname(b)];
        if (count.mate1_count > 1 || count.mate2_count > 1) {
            continue;
        }

        while(i < intervals.size() && (b->core.tid > intervals[i].tid ||
              (b->core.tid == intervals[i].tid &&
               start_pos > intervals[i].locus.max_end))) {
            process_locus(counts, intervals[i].locus, rs, stranded);
            rs.clear();
            ++i;
        }
        if (i >= intervals.size()) break;

        if (b->core.tid < intervals[i].tid) continue;
        if (b->core.tid == intervals[i].tid) {
            if (start_pos < intervals[i].locus.min_start) continue;

            if (start_pos <= intervals[i].locus.max_end &&
                end_pos >= intervals[i].locus.min_start) {
                rs.add_alignment(idx, b);
                continue;
            }
        }
    }

    samclose(sam_file);
    bam_destroy1(b);

    // print output
    fprintf(output_file, "id\tcount\n");
    for (TrieMapIterator<unsigned long> c(counts);
         c != TrieMapIterator<unsigned long>(); ++c) {
        fprintf(output_file, "%s\t%lu\n", c->first, *c->second);
    }

    if (output_file != stdout) fclose(output_file);

    return EXIT_SUCCESS;
}



