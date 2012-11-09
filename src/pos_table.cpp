
#include "pos_table.hpp"

PosTable::PosTable()
{
}


PosTable::~PosTable()
{
}


size_t PosTable::size()
{
    size_t n = 0;
    std::vector<PosSubtable>::iterator subtable;
    for (subtable = subtables.begin(); subtable != subtables.end(); ++subtable) {
        int strand;
        for (strand = 0; strand < 2; ++strand) {
            std::vector<PosTableVal>::iterator i;
            for (i = subtable->values[strand].begin();
                 i != subtable->values[strand].end();
                 ++i) {
                ++n;
            }
        }
    }

    return n;
}


void PosTable::add(bam1_t* b, samfile_t* bam_f)
{
    if (subtables.empty() || subtables.back().tid != b->core.tid) {
        subtables.push_back(PosSubtable());
        subtables.back().tid = b->core.tid;
        subtables.back().seqname = bam_f->header->target_name[b->core.tid];
    }

    PosSubtable& subtable = subtables.back();

    int strand = bam1_strand(b);
    int32_t pos;
    if (strand) pos = bam_calend(&b->core, bam1_cigar(b)) - 1;
    else        pos = b->core.pos;

    if (subtable.values[strand].empty() || subtable.values[strand].back().pos != pos) {
        subtable.values[strand].push_back(PosTableVal());
        subtable.values[strand].back().pos = pos;
        subtable.values[strand].back().count = 1;
    }
    else {
        subtable.values[strand].back().count++;
    }
}


void PosTable::dump(std::vector<ReadPos>& positions)
{
    size_t k = 0;
    std::vector<PosSubtable>::iterator subtable;
    for (subtable = subtables.begin(); subtable != subtables.end(); ++subtable) {
        int strand;
        for (strand = 0; strand < 2; ++strand) {
            std::vector<PosTableVal>::iterator i;
            for (i = subtable->values[strand].begin();
                 i != subtable->values[strand].end();
                 ++i) {
                if (k < positions.size()) {
                    positions[k].seqname = subtable->seqname;
                    positions[k].strand  = strand;
                    positions[k].pos     = i->pos;
                    positions[k].count   = i->count;
                    ++k;
                }
                else {
                    goto finish_pos_table_dump;
                }
            }
        }
    }
finish_pos_table_dump:
    return;
}


