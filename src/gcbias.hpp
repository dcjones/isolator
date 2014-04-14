
#ifndef ISOLATOR_GC_BIAS
#define ISOLATOR_GC_BIAS

#include <vector>

#include "common.hpp"
#include "pos_table.hpp"
#include "samtools/faidx.h"
#include "intervals.hpp"
#include "seqbias/sequencing_bias.hpp"


class GCBias
{
public:
	GCBias(const char* ref_fn, PosTable& T,
		   pos_t median_frag_len,
           sequencing_bias* seqbias[2]);
	~GCBias();
	double get_bias(double gc);

private:
	std::vector<double> bins;
	std::vector<double> bin_bias;
};

#endif
