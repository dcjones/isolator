
#ifndef ISOLATOR_SUMMARIZE_HPP
#define ISOLATOR_SUMMARIZE_HPP

// Read and produce plain-text tables from the HDF5 output generated
// by "analyze".

#include <boost/numeric/ublas/matrix.hpp>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>


class Summarize
{
	public:
		Summarize(const char* filename);
		~Summarize();

		// pre-baked summarizations
		void median_condition_tgroup_expression(FILE* output);
		void median_experiment_tgroup_sd(FILE* output);

        void median_transcript_expression(FILE* output);
        void median_gene_expression(FILE* output);


	private:
        void median_transcript_expression(
                boost::numeric::ublas::matrix<float>& Q);

		hid_t h5_file;

        // number of samples
        size_t K;

		// number of transcripts
		size_t N;

		// number of tgroups
		size_t T;

		// numebr of conditions
		size_t C;

		// transcript_id indexed by tid
		std::vector<std::string> transcript_ids;

		// gene_id indexed by tid
		std::vector<std::string> gene_ids;

		// tgroup indeed by tid
		std::vector<unsigned int> tgroup;
};


#endif
