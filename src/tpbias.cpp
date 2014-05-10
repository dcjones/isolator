
#include <cmath>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

#include "fastmath.hpp"
#include "nlopt/nlopt.h"
#include "tpbias.hpp"


typedef std::pair<pos_t, pos_t> TPDistPair;
typedef boost::unordered_map<pos_t, double> DistNormMap;


static double geometric_cdf_upper(double p, unsigned int k)
{
    return std::pow(1 - p, k);
}


static double geometric_logcdf_upper(double p, unsigned int k)
{
    return k * fastlog(1 - p);
}


static double tpbias_loglikelihood(double p,
                                   pos_t maxtlen,
                                   DistNormMap& distnorm,
                                   const std::vector<TPDistPair>& xs)
{
    double norm = 0.0;
    double ppower = 1.0;
    for (pos_t tlen = 0; tlen <= maxtlen; ++tlen) {
        ppower *= 1 - p;
        norm += ppower;

        DistNormMap::iterator i = distnorm.find(tlen);
        if (i != distnorm.end()) {
            i->second = fastlog(norm);
        }
    }

    double ll = 0.0;
    BOOST_FOREACH (const TPDistPair& x, xs)  {
        ll += geometric_logcdf_upper(p, x.first) - distnorm[x.second];
    }

    return ll;
}


struct TPBiasTrainData
{
    TPBiasTrainData(pos_t maxtlen,
                    DistNormMap& distnorm,
                    const std::vector<TPDistPair>& xs)
        : maxtlen(maxtlen)
        , distnorm(distnorm)
        , xs(xs)
    {
    }

    pos_t maxtlen;
    DistNormMap& distnorm;
    const std::vector<TPDistPair>& xs;
};


static double tpbias_nlopt_objective(unsigned int _n, const double* _x,
                                     double* _grad, void* data)
{
    UNUSED(_n);
    UNUSED(_grad);

    TPBiasTrainData* traindata =
        reinterpret_cast<TPBiasTrainData*>(data);

    return tpbias_loglikelihood(_x[0], traindata->maxtlen,
                                traindata->distnorm, traindata->xs);
}


TPBias::TPBias(const std::vector<std::pair<pos_t, pos_t> >& tpdists)
{
    // normalizing constants for fragment disributions given tlen
    boost::unordered_map<pos_t, double> distnorm;

    pos_t maxtlen = 0;
    BOOST_FOREACH (const TPDistPair& tpdist_tlen, tpdists) {
        distnorm[tpdist_tlen.second] = 0.0;
        maxtlen = std::max<pos_t>(maxtlen, tpdist_tlen.second);
    }

    TPBiasTrainData traindata(maxtlen, distnorm, tpdists);

    nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX, 1);
    double lower_limit = 0.0;
    double upper_limit = 1e-4;
    nlopt_set_lower_bounds(opt, &lower_limit);
    nlopt_set_upper_bounds(opt, &upper_limit);
    nlopt_set_max_objective(opt, tpbias_nlopt_objective,
                            reinterpret_cast<void*>(&traindata));
    nlopt_set_ftol_abs(opt, 1e-4);

    //double xtol_abs = 1e-9;
    //nlopt_set_xtol_abs(opt, &xtol_abs);

    p = 1e-5;
    double maxf;
    nlopt_result result = nlopt_optimize(opt, &p, &maxf);
    if (result < 0) {
        Logger::warn("Failed to fit 3' bias model.");
        p = 0.0;
    }

    nlopt_destroy(opt);
}


double TPBias::get_bias(pos_t k)
{
    return geometric_cdf_upper(p, k);
}


