
#include "intervals.hpp"

Interval::Interval()
    : start(-1)
    , end(-1)
    , strand(strand_na)
{
}


Interval::Interval(const char* seqname, pos_t start, pos_t end, strand_t strand)
    : seqname(seqname)
    , start(start)
    , end(end)
    , strand(strand)
{
}


Interval::Interval(const Interval& other)
    : seqname(other.seqname)
    , start(other.start)
    , end(other.end)
    , strand(other.strand)
{
}


pos_t Interval::length() const
{
    return end - start + 1;
}


bool Interval::operator < (const Interval& other) const
{
    if (seqname != other.seqname) return seqname < other.seqname;
    if (strand  != other.strand)  return strand < other.strand;
    if (start   != other.start)   return start < other.start;
    if (end     != other.end)     return end < other.end;
    return false;
}


bool Interval::operator == (const Interval& other) const
{
    return seqname == other.seqname &&
           strand  == other.strand &&
           start   == other.start &&
           end     == other.end;
}
