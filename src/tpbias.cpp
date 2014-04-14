
#include <cmath>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

#include "tpbias.hpp"
#include "fastmath.hpp"


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


TPBias::TPBias(const std::vector<std::pair<pos_t, pos_t> >& tpdists)
{
    // normalizing constants for fragment disributions given tlen
    boost::unordered_map<pos_t, double> distnorm;

    pos_t maxtlen = 0;
    BOOST_FOREACH (const TPDistPair& tpdist_tlen, tpdists) {
        distnorm[tpdist_tlen.second] = 0.0;
        maxtlen = std::max<pos_t>(maxtlen, tpdist_tlen.second);
    }

    // binary search of the maximum likelihood p
    double p0 = 0.01;
    double step = 2.0;
    double ll0 = tpbias_loglikelihood(p0, maxtlen, distnorm, tpdists);

    // Janky hillclimbing. Assume max_p is suboptimal and ll should increase
    // monotonically as we move towards 0.0, and the optimal answer. So we
    // just keep dividing by two until the likelihood doesn't increase.

    while (step > 1.1) {
        double p = p0 / step;
        double ll = tpbias_loglikelihood(p, maxtlen, distnorm, tpdists);
        if (ll > ll0) {
            p0 = p;
            ll0 = ll;
        }
        else {
            step *= 0.9;
        }
    }

    this->p = p0;
}


double TPBias::get_bias(pos_t k)
{
    return geometric_cdf_upper(p, k);
}


