
#ifndef ISOLATOR_GC_HPP
#define ISOLATOR_GC_HPP

#include <vector>

#include "transcripts.hpp"
#include "loess/loess.h"


class GCCorrection
{
    public:
        GCCorrection(TranscriptSet& ts, const float* transcript_gc);
        ~GCCorrection();

        void correct(double* expr);

    private:
        TranscriptSet& ts;

        const float* transcript_gc;

        // number of entrie (tgroups) in xs, ys
        size_t n;

        // input: weighted gc content
        std::vector<double> xs;

        // response: log tgroup expression
        std::vector<double> ys;

        // transcript gc
        std::vector<double> tgc;

        // loess predictions
        std::vector<double> fit, se_fit;

        loess_struct lo;
        pred_struct pre;
};



#endif

