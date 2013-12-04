
#ifndef ISOLATOR_INTERVALS_HPP
#define ISOLATOR_INTERVALS_HPP

#include "common.hpp"

struct Interval
{
    Interval();
    Interval(const char* seqname, pos_t start, pos_t end, strand_t strand);
    Interval(SeqName seqname, pos_t start, pos_t end, strand_t strand);
    Interval(const Interval&);

    bool operator < (const Interval&) const;
    bool operator == (const Interval&) const;

    pos_t length() const;

    SeqName seqname;
    pos_t start, end;
    strand_t strand;
};


struct IntervalPair
{
    IntervalPair();
    IntervalPair(const Interval& first, const Interval& second);

    bool operator < (const IntervalPair&) const;
    bool operator == (const IntervalPair&) const;

    Interval first, second;
};


#endif

