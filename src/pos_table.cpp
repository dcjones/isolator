
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


void PosTable::add(int32_t tid, pos_t pos, int strand,
                   pos_t start, pos_t end, samfile_t* bam_f)
{
    if (subtables.empty() || subtables.back().tid != tid) {
        subtables.push_back(PosSubtable());
        subtables.back().tid = tid;
        subtables.back().seqname = bam_f->header->target_name[tid];
    }

    PosSubtable& subtable = subtables.back();

    if (subtable.values[strand].empty() ||
        subtable.values[strand].back().pos != pos) {
        subtable.values[strand].push_back(PosTableVal());
        subtable.values[strand].back().pos = pos;
        subtable.values[strand].back().start = start;
        subtable.values[strand].back().end = end;
        subtable.values[strand].back().count = 1;
    }
    else {
        subtable.values[strand].back().count++;
    }
}


void PosTable::dump(std::vector<ReadPos>& positions, size_t max_size)
{
    std::vector<PosSubtable>::iterator subtable;
    for (subtable = subtables.begin(); subtable != subtables.end(); ++subtable) {
        int strand;
        for (strand = 0; strand < 2; ++strand) {
            std::vector<PosTableVal>::iterator i;
            for (i = subtable->values[strand].begin();
                 i != subtable->values[strand].end();
                 ++i) {
                if (positions.size() < max_size) {
                    positions.push_back(
                            ReadPos(subtable->seqname, strand, i->pos,
                                    i->start, i->end, i->count));
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


