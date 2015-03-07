
#include <algorithm>
#include <boost/foreach.hpp>

#include "constants.hpp"
#include "fragment_model.hpp"
#include "fragbias.hpp"


FragBias::FragBias(FragmentModel* fm)
{
    // lengths used for interpolation
    const pos_t lengths[] = {
        25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300,
        325, 350, 375, 400, 425, 450, 475, 500, 550, 600, 650, 700,
        750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1500, 1600,
        1700, 1800, 1900, 2000, 2500, 3000, 3500, 4000, 4500, 5000
    };

    this->lengths.resize(sizeof(lengths) / sizeof(pos_t));
    std::copy(lengths, lengths + this->lengths.size(), this->lengths.begin());

    // assumed fragmentation probability
    this->fm = fm;
    p = (1.0 / fm->frag_len_med()) / 2.0;
    double omp = 1.0 - p;
    max_non_frag_run = std::max<pos_t>(1, lround(log(1e-5) / log(omp)));

    non_frag_run_pr_sum.resize(max_non_frag_run);
    non_frag_run_pr.resize(max_non_frag_run);

    non_frag_run_pr[0] = 1.0;
    non_frag_run_pr_sum[0] = 1.0;
    for (pos_t m = 1; m < max_non_frag_run; ++m) {
        non_frag_run_pr_sum[m] = non_frag_run_pr_sum[m - 1] + pow(omp, m - 1) / (double) m;
        non_frag_run_pr[m] = pow(omp, m - 1);
    }

    BOOST_FOREACH (pos_t l, this->lengths) {
        std::vector<float> fb(l);
        compute(&fb.at(0), 0, l, l);
        bias.push_back(fb);
    }
}


void FragBias::compute(float* out, pos_t from, pos_t to, pos_t l)
{
        // Four cases:
        for (pos_t u = from; u < to; ++u) out[u] = 0.0;

#if 0
        // enrichment of end positions assuming the end position is
        // primed first
        // Case 1: no fragmentation
        double c = (1.0 / l) * (l - 1 < max_non_frag_run ? non_frag_run_pr[l - 1] : 0.0);
        for (pos_t v = from; v < to; ++v) {
            double select_pr = fm->frag_len_c(v) / v;
            out[v] = c * select_pr;
        }

        // Case 2: read contained in a fragment at the left end of the
        //         transcript
        for (pos_t v = from; v < to; ++v) {
            for (pos_t j = v; j < std::min<pos_t>(l, v + max_non_frag_run); ++j) {
                double select_pr = fm->frag_len_c(v) / v;
                double run_pr = j - 1 < max_non_frag_run ? non_frag_run_pr[j - 1] : 0.0;
                out[v] += p * run_pr * (1.0 / j) * select_pr;
            }
        }

        // Case 3: read contained a fragment at the right end
        for (pos_t v = from; v < to; ++v) {
            for (pos_t i = std::max<pos_t>(0, v - max_non_frag_run); i < v; ++i) {
                double select_pr = fm->frag_len_c(v - i) / (v - i);
                double run_pr = l - i - 1 < max_non_frag_run ? non_frag_run_pr[l - i - 1] : 0.0;
                out[v] += p * run_pr * (1.0 / ((l - i))) * select_pr;
            }
        }

        // Case 4: read flanked by fragmentation breakpoints
        for (pos_t v = from; v < to; ++v) {
            for (pos_t h = 1; h < std::min<pos_t>(max_non_frag_run, v); ++h) {
                double run_pr1 = non_frag_run_pr_sum[
                    std::min<pos_t>((l - v)+ h, max_non_frag_run - 1)];
                double run_pr2 = non_frag_run_pr_sum[
                    std::min<pos_t>(h, max_non_frag_run - 1)];

                out[v] += p * p * (fm->frag_len_c(h) / h) *
                    (run_pr1 - run_pr2);
            }
        }
#endif

        // enrichment of start positions assuming the start position
        // is primed first.
        // Case 1: no fragmentation
        double c = (1.0 / l) * (l - 1 < max_non_frag_run ? non_frag_run_pr[l - 1] : 0.0);
        for (pos_t u = from; u < to; ++u) {
            double select_pr = fm->frag_len_c(l - u + 1) / (l - u + 1);
            out[u] = c * select_pr;
        }

        // Case 2: read contained in a fragment at the left end of the
        //         transcript
        for (pos_t u = from; u < to; ++u) {
            for (pos_t j = u + 1; j < std::min<pos_t>(l, u + max_non_frag_run); ++j) {
                double select_pr = fm->frag_len_c(j - u) / (j - u);
                double run_pr = j - 1 < max_non_frag_run ? non_frag_run_pr[j - 1] : 0.0;
                out[u] += p * run_pr * (1.0 / j) * select_pr;
            }
        }

        // Case 3: read contained a fragment at the right end
        for (pos_t u = from; u < to; ++u) {
            for (pos_t i = std::max<pos_t>(0, u - max_non_frag_run); i < u; ++i) {
                double select_pr = fm->frag_len_c(l - u + 1) / (l - u + 1);
                double run_pr = l - i - 1 < max_non_frag_run ? non_frag_run_pr[l - i - 1] : 0.0;
                out[u] += p * run_pr * (1.0 / ((l - i))) * select_pr;
            }
        }

        // Case 4: read flanked by fragmentation breakpoints
        for (pos_t u = from; u < to; ++u) {
            for (pos_t h = 1; h < std::min<pos_t>(max_non_frag_run, l - u); ++h) {
                double run_pr1 = non_frag_run_pr_sum[
                    std::min<pos_t>(u + h, max_non_frag_run - 1)];
                double run_pr2 = non_frag_run_pr_sum[
                    std::min<pos_t>(h, max_non_frag_run - 1)];

                out[u] += p * p * (fm->frag_len_c(h) / h) *
                    (run_pr1 - run_pr2);
            }
        }

        for (pos_t u = from; u < to; ++u) {
            out[u] *= constants::fragbias_scale;
        }
}


void FragBias::get_bias(float* out, pos_t l)
{
    unsigned int a;
    // extrapolate downwards
    if (l < lengths.front()) {
        a = 0;
    }
    // extrapolate upwards
    else if (l >= lengths.back()) {
        a = lengths.size() - 2;
    }
    // interpolate
    else {
        std::vector<pos_t>::iterator a_it =
            std::upper_bound(lengths.begin(), lengths.end(), l);
        a = a_it - lengths.begin() - 1;
    }

    unsigned int b = a + 1;

    pos_t al = lengths[a];
    pos_t bl = lengths[b];

    double z = (double) (l - al) / (double) (bl - al);

    // extrapolation can be very very wrong, so we just explicitly compute the
    // ends and copy the middles.
    pos_t endlen = 0;
    if (b == lengths.size() - 1) {
        z = 0.0;
        endlen = constants::fragbias_endlen;
    }
    pos_t to = std::min<pos_t>(l, endlen);
    compute(out, 0, to, l);
    compute(out, std::max<pos_t>(to, l - endlen), l, l);

    for (pos_t u = endlen; u < l - endlen; ++u) {
        double r = (double) u / (double) (l - 1);
        pos_t au = lround(r * (al - 1));
        pos_t bu = lround(r * (bl - 1));
        double slope = bias[b][bu] - bias[a][au];

        out[u] = bias[a][au] + z * slope;
    }
}


