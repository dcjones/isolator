
#ifndef ISOLATOR_GC_HPP
#define ISOLATOR_GC_HPP

#include <vector>

#include "transcripts.hpp"
#include "loess/loess.h"


class GCCorrection
{
    public:
        GCCorrection(TranscriptSet& ts, const double* transcript_gc);
        ~GCCorrection();

        void correct(double* expr);
        void adjustments(const double* expr, double* weights);

    private:
        TranscriptSet& ts;

        const double* transcript_gc;

        // number of entrie (tgroups) in xs, ys
        size_t n;

        // work vectors for computing gene expression and weighted gc content
        std::vector<double> gene_expr, gene_gc;

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

