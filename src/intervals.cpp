
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

