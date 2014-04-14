
#ifndef ISOLATOR_TPBIAS
#define ISOLATOR_TPBIAS

#include <vector>

#include "common.hpp"

class TPBias
{
    public:
        TPBias(const std::vector<std::pair<pos_t, pos_t> >& tpdists);
        double get_bias(pos_t k);

        double p;
};


#endif
