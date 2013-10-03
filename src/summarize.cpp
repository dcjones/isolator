
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "summarize.hpp"
#include "logger.hpp"


using namespace boost::numeric::ublas;


Summarize::Summarize(const char* filename)
{
	h5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (h5_file < 0) {
		Logger::abort("Failed to open HDF5 file %s", filename);
	}

	// read transcript information
	hid_t dataset;

	dataset = H5Dopen2(h5_file, "/transcript_id", H5P_DEFAULT);
	if (dataset < 0) {
		Logger::abort("Could not oepen the transcript_id dataset.");
	}

	hid_t dataspace = H5Dget_space(dataset);

	hsize_t dims[1];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	N = dims[0];
	Logger::info("%lu transcripts", (unsigned long) N);

	hid_t datatype = H5Tcopy(H5T_C_S1);
	H5Tset_size(datatype, H5T_VARIABLE);

	const char** string_data = new const char* [N];

	// read transcript_id dataset
	herr_t status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL,
		                    H5P_DEFAULT, string_data);
	if (status < 0) {
		Logger::abort("Failed to read the transcript_id dataset.");
	}

	transcript_ids.reserve(N);
	for (size_t i = 0; i < N; ++i) {
		transcript_ids.push_back(std::string(string_data[i]));
	}

	H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, string_data);
	H5Dclose(dataset);

	// read gene_id dataset
	dataset = H5Dopen2(h5_file, "/gene_id", H5P_DEFAULT);
	if (dataset < 0) {
		Logger::abort("Could not open the gene_id dataset.");
	}

	status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL,
		             H5P_DEFAULT, string_data);
	if (status < 0) {
		Logger::abort("Failed to read the gene_id dataset.");
	}

	gene_ids.reserve(N);
	for (size_t i = 0; i < N; ++i) {
		gene_ids.push_back(std::string(string_data[i]));
	}

	H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, string_data);
	H5Dclose(dataset);

	delete [] string_data;
	H5Tclose(datatype);

	// read tgroups
	dataset = H5Dopen2(h5_file, "/tgroup", H5P_DEFAULT);
	if (dataset < 0) {
		Logger::abort("Could not open the tgroup dataset.");
	}

	tgroup.resize(N);
	status = H5Dread(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
		             H5P_DEFAULT, &tgroup.at(0));
	if (status < 0) {
		Logger::abort("Failed to read the tgroup dataset");
	}

	T = *std::max_element(tgroup.begin(), tgroup.end()) + 1;
	Logger::info("%lu tgroups", (unsigned long) T);

	H5Dclose(dataset);
    H5Sclose(dataspace);

    // figure out K (number of samples)
    dataset = H5Dopen2(h5_file, "/transcript_quantification", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /transcript_quant dataset");
    }

    dataspace = H5Dget_space(dataset);
    hsize_t dims3[3]; // dims are num_samples, K, N
    H5Sget_simple_extent_dims(dataspace, dims3, NULL);
    K = dims3[1];

    H5Dclose(dataset);
    H5Sclose(dataspace);

}


void Summarize::median_condition_tgroup_expression(FILE* output)
{
	typedef std::vector<std::set<std::string> > string_set_vector;

	string_set_vector tgroup_gene_ids(T);
	string_set_vector tgroup_transcript_ids(T);

	for (size_t i = 0; i < N; ++i) {
		tgroup_gene_ids[tgroup[i]].insert(gene_ids[i]);
		tgroup_transcript_ids[tgroup[i]].insert(transcript_ids[i]);
	}

	hid_t dataset = H5Dopen2(h5_file, "/condition/tgroup_mean", H5P_DEFAULT);
	if (dataset < 0) {
		Logger::abort("Failed to open the /condition/tgroup_mean dataset");
	}

	hid_t dataspace = H5Dget_space(dataset);
	hsize_t dims[3]; // dims are: num_samples, C, T
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	C = dims[1];
	size_t num_samples = dims[0];

	hsize_t mem_dataspace_dims[1] = { T };
	hsize_t mem_dataspace_start[1] = { 0 };
	hid_t mem_dataspace = H5Screate_simple(1, mem_dataspace_dims, NULL);
	H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, mem_dataspace_start,
		                NULL, mem_dataspace_dims, NULL);

	// temporary vector for computing medians
	std::vector<float> median_work(num_samples);

	hsize_t file_dataspace_dims[3] = {1, 1, T};
	hsize_t file_dataspace_start[3] = {0, 0, 0};

	// tgroup_mean posterior medians indexed by tgroup and condition
	matrix<float> tgroup_medians(T, C);

	// all tgroup_mean data indexed by sample_num and tgroup
	std::vector<std::vector<float> > condition_data(num_samples);
	for (size_t i = 0; i < num_samples; ++i) {
		condition_data[i].resize(T);
	}

	for (size_t i = 0; i < C; ++i) {
		for (size_t j = 0; j < num_samples; ++j) {
			file_dataspace_start[0] = j;
			file_dataspace_start[1] = i;
			H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
				                file_dataspace_start, NULL,
				                file_dataspace_dims, NULL);

			herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT,
				                    mem_dataspace, dataspace, H5P_DEFAULT,
				                    &condition_data[j].at(0));

			if (status < 0) {
				Logger::abort("Reading /condition/tgroup_mean dataset failed.");
			}
		}

		for (size_t j = 0; j < T; ++j) {
			for (size_t k = 0; k < num_samples; ++k) {
				median_work[k] = condition_data[k][j];
			}
			std::sort(median_work.begin(), median_work.end());
			tgroup_medians(j, i) = median_work[median_work.size() / 2];
		}
	}

	H5Sclose(mem_dataspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);

	fprintf(output, "transcript_ids\tgene_ids");
	for (unsigned int i = 0; i < C; ++i) {
		fprintf(output, "\tcondition%u_mean_expr", i + 1);
	}
	fputc('\n', output);

	for (size_t i = 0; i < T; ++i) {
		bool firstid = true;
		BOOST_FOREACH (const std::string& transcript_id, tgroup_transcript_ids[i]) {
			fprintf(output, firstid ? "%s" : ",%s", transcript_id.c_str());
			firstid = false;
		}
		fputc('\t', output);

		firstid = true;
		BOOST_FOREACH (const std::string& gene_id, tgroup_gene_ids[i]) {
			fprintf(output, firstid ? "%s" : ",%s", gene_id.c_str());
			firstid = false;
		}

		for (size_t j = 0; j < C; ++j) {
			fprintf(output, "\t%f", tgroup_medians(i, j));
		}
		fputc('\n', output);
	}
}


void Summarize::median_experiment_tgroup_sd(FILE* output)
{
	typedef std::vector<std::set<std::string> > string_set_vector;

	string_set_vector tgroup_gene_ids(T);
	string_set_vector tgroup_transcript_ids(T);

	for (size_t i = 0; i < N; ++i) {
		tgroup_gene_ids[tgroup[i]].insert(gene_ids[i]);
		tgroup_transcript_ids[tgroup[i]].insert(transcript_ids[i]);
	}

	hid_t dataset = H5Dopen2(h5_file, "/experiment/tgroup_sd", H5P_DEFAULT);
	if (dataset < 0) {
		Logger::abort("Failed to open the /experiment/tgroup_sd dataset");
	}

	hid_t dataspace = H5Dget_space(dataset);
	hsize_t dims[2]; // dims are: num_samples, T
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	size_t num_samples = dims[0];
	T = dims[1];

	hsize_t mem_dataspace_dims[1] = {T};
	hsize_t mem_dataspace_start[1] = {0};
	hid_t mem_dataspace = H5Screate_simple(1, mem_dataspace_dims, NULL);
	H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, mem_dataspace_start,
		                NULL, mem_dataspace_dims, NULL);

	std::vector<std::vector<float> > data(num_samples);
	for (size_t i = 0; i < num_samples; ++i) {
		data[i].resize(T);
	}

	hsize_t file_dataspace_start[2] = {0, 0};
	hsize_t file_dataspace_dims[2] = {1, T};

	for (size_t i = 0; i < num_samples; ++i) {
		file_dataspace_start[0] = i;
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			                file_dataspace_start, NULL,
			                file_dataspace_dims, NULL);

		herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT, mem_dataspace,
			                    dataspace, H5P_DEFAULT, &data[i].at(0));

		if (status < 0) {
			Logger::abort("Reading the /experiment/tgroup_sd dataset failed.");
		}
	}

	H5Sclose(mem_dataspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);

	// temporary vector for computing medians
	std::vector<float> median_work(num_samples);

	fprintf(output, "transcript_ids\tgene_ids\tsd\n");
	for (size_t i = 0; i < T; ++i) {
		bool firstid = true;
		BOOST_FOREACH (const std::string& transcript_id, tgroup_transcript_ids[i]) {
			fprintf(output, firstid ? "%s" : ",%s", transcript_id.c_str());
			firstid = false;
		}
		fputc('\t', output);

		firstid = true;
		BOOST_FOREACH (const std::string& gene_id, tgroup_gene_ids[i]) {
			fprintf(output, firstid ? "%s" : ",%s", gene_id.c_str());
			firstid = false;
		}

		for (size_t j = 0; j < num_samples; ++j) {
			median_work[j] = data[j][i];
		}
		std::sort(median_work.begin(), median_work.end());
		fprintf(output, "\t%f\n", median_work[median_work.size() / 2]);
	}
}


void Summarize::median_transcript_expression(matrix<float>& Q)
{
    hid_t dataset = H5Dopen2(h5_file, "/transcript_quantification", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /transcript_quant dataset");
    }

    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[3]; // dims are num_samples, K, N
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t num_samples = dims[0];
    size_t K = dims[1];

    hsize_t file_dataspace_start[3] = {0, 0, 0};
    hsize_t file_dataspace_dims[3] = {num_samples, 1, N};

    hsize_t mem_dataspace_dims[2] = {num_samples, N};
    hsize_t mem_dataspace_start[2] = {0, 0};
    hid_t mem_dataspace = H5Screate_simple(2, mem_dataspace_dims, NULL);
    H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, mem_dataspace_start,
                        NULL, mem_dataspace_dims, NULL);

    matrix<float> Qi(num_samples, N);

    for (size_t i = 0; i < K; ++i) {
        herr_t status;
        file_dataspace_start[1] = i;

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            file_dataspace_start, NULL,
                            file_dataspace_dims, NULL);

        status = H5Dread(dataset, H5T_NATIVE_FLOAT, mem_dataspace,
                         dataspace, H5P_DEFAULT, &Qi.data()[0]);

        if (status < 0) {
            Logger::abort("Reading the /transcript_quantification dataset failed.");
        }

        for (size_t j = 0; j < N; ++j) {
            matrix_column<matrix<float> > col(Qi, j);
            std::sort(col.begin(), col.end());
            float med = col[col.size() / 2];
            Q(i, j) = med;
        }
    }

	H5Sclose(mem_dataspace);
	H5Sclose(dataspace);
	H5Dclose(dataset);
}


void Summarize::median_transcript_expression(FILE* output)
{
    matrix<float> Q(K, N);
    median_transcript_expression(Q);

    fprintf(output, "transcript_id\tgene_id");
    for (unsigned int i = 0; i < K; ++i) {
        fprintf(output, "\tsample_%u", i);
    }
    fputc('\n', output);

    for (size_t i = 0; i < N; ++i) {
        fprintf(output, "%s\t%s",
                transcript_ids[i].c_str(),
                gene_ids[i].c_str());
        for (size_t j = 0; j < K; ++j) {
            fprintf(output, "\t%e", Q(j, i));
        }
        fputc('\n', output);
    }
}


void Summarize::median_gene_expression(FILE* output)
{
    matrix<float> Q(K, N);
    median_transcript_expression(Q);

    std::map<std::string, std::vector<size_t> > gid_to_tids;
    for (size_t i = 0; i < N; ++i) {
        gid_to_tids[gene_ids[i]].push_back(i);
    }

    fprintf(output, "gene_id");
    for (unsigned int i = 0; i < K; ++i) {
        fprintf(output, "\tsample_%u", i);
    }
    fputc('\n', output);

    boost::numeric::ublas::vector<float> margin(K);

    typedef std::pair<std::string, std::vector<size_t> > item_t;
    BOOST_FOREACH (const item_t& item, gid_to_tids) {
        fprintf(output, "%s", item.first.c_str());
        std::fill(margin.begin(), margin.end(), 0.0);
        BOOST_FOREACH (size_t tid, item.second) {
            margin += matrix_column<matrix<float> >(Q, tid);
        }

        for (size_t j = 0; j < K; ++j) {
            fprintf(output, "\t%e", margin[j]);
        }
        fputc('\n', output);
    }
}


Summarize::~Summarize()
{
	H5Fclose(h5_file);
}

