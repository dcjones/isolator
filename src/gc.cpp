
#include <boost/foreach.hpp>

#include "constants.hpp"
#include "gc.hpp"


GCCorrection::GCCorrection(TranscriptSet& ts, const float* transcript_gc)
    : ts(ts)
    , transcript_gc(transcript_gc)
{
    n = 0;
    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        n = std::max<size_t>(n, t->tgroup + 1);
    }

    xs.resize(n);
    ys.resize(n);
    fit.resize(ts.size());
    se_fit.resize(ts.size());
    tgc.resize(ts.size());

    loess_setup2(&xs.at(0), &ys.at(0), n, 1, &lo);
    lo.model.span = constants::gc_loess_smoothing;
    lo.model.family = "symmetric";
    lo.model.degree = 2;
}


GCCorrection::~GCCorrection()
{
    loess_free_mem2(&lo);
}


void GCCorrection::correct(double* expr)
{
    std::fill(xs.begin(), xs.end(), 0.0);
    std::fill(xs.begin(), xs.end(), 0.0);
    for (TranscriptSet::iterator t = ts.begin(); t != ts.end(); ++t) {
        if (expr[t->id] == 0.0 ||
            t->exonic_length() < 200 ||
            transcript_gc[t->id] < 0.20 ||
            transcript_gc[t->id] > 0.70) {
            continue;
        }

        xs[t->id] += expr[t->id] * transcript_gc[t->id];
        ys[t->id] += expr[t->id];
    }

    double max_gc = 0.0, min_gc = 1.0;
    for (size_t i = 0; i < n; ++i) {
        if (ys[i] > 0.0) {
            xs[i] /= ys[i];
            ys[i] = log(ys[i]);
        }
        else {
            ys[i] = log(constants::zero_eps);
        }
        max_gc = std::max<double>(max_gc, xs[i]);
        min_gc = std::min<double>(min_gc, xs[i]);
    }

    // normalize to a arbitary point in the loess curve
    double ref_gc = 0.5;
    double ref, ref_se;
    predict2(&lo, &ref_gc, &ref, &ref_se, 1, FALSE);

    // normalize samples
    std::copy(transcript_gc, transcript_gc + ts.size(), tgc.begin());
    for (size_t i = 0; i < ts.size(); ++i) {
        tgc[i] = std::min(max_gc, std::max(min_gc, tgc[i]));
    }

    predict2(&lo, &tgc.at(0), &fit.at(0), &se_fit.at(0), ts.size(), FALSE);

    double z = 0.0;
    for (size_t i = 0; i < ts.size(); ++i) {
        expr[i]  = pow(expr[i], ref / fit[i]);
        z += expr[i];
    }

    for (size_t i = 0; i < ts.size(); ++i) {
        expr[i] /= z;
    }
}


