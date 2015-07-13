
#include "pos_table.hpp"


PosTable::PosTable()
{}


PosTable::~PosTable()
{}


bool PosSubtable::inc(pos_t pos, pos_t start, pos_t end, int strand)
{
    int s = strand == 0 ? 0 : 1;

    boost::unordered_map<pos_t, PosTableVal>::iterator i = table[s].find(pos);
    if (i == table[s].end()) {
        table[s].insert(std::pair<pos_t, PosTableVal>(
            pos, PosTableVal(pos, start, end, 1)));
        return true;
    }
    else {
        i->second.count++;
        return false;
    }
}


size_t PosTable::size()
{
    size_t n = 0;
    std::vector<PosSubtable>::iterator subtable;
    for (subtable = subtables.begin(); subtable != subtables.end(); ++subtable) {
        int strand;
        for (strand = 0; strand < 2; ++strand) {
            n += subtable->table[strand].size();
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
    subtable.inc(pos, start, end, strand);
}


void PosTable::dump(std::vector<ReadPos>& positions, size_t max_size)
{
    std::vector<PosSubtable>::iterator subtable;
    for (subtable = subtables.begin(); subtable != subtables.end(); ++subtable) {
        int strand;
        for (strand = 0; strand < 2; ++strand) {
            boost::unordered_map<pos_t, PosTableVal>::iterator i;
            for (i = subtable->table[strand].begin();
                 i != subtable->table[strand].end();
                 ++i) {
                if (positions.size() < max_size) {
                    positions.push_back(
                            ReadPos(subtable->seqname, strand, i->second.pos,
                                    i->second.start, i->second.end, i->second.count));
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


