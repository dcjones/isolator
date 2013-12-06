
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <intervals.hpp>
#include <set>

#include "summarize.hpp"
#include "logger.hpp"


using namespace boost::numeric::ublas;


// Compute log(exp(x) + exp(y)), avoiding overflow/underflow.
static double logaddexp(double x, double y)
{
    double u = x - y;
    if (u > 0.0) {
        return x + log1p(exp(-u));
    }else if (u <= 0.0) {
        return y + log1p(exp(u));
    }else  {
        return x + y;
    }
}


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

    // construct tgroup_tids
    tgroup_tids.resize(T);
    for (size_t i = 0; i < N; ++i) {
        tgroup_tids[tgroup[i]].push_back(i);
    }

    // construct spliced_tgroup_indexes
    for (size_t i = 0; i < T; ++i) {
        if (tgroup_tids[i].size() > 1) {
            spliced_tgroup_indexes.push_back(i);
        }
    }

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
    size_t C = dims[1];
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

void Summarize::tgroup_fold_change(FILE* output,
                                   unsigned int condition_a,
                                   unsigned int condition_b)
{
    // TODO: pass these in
    double upper_quantile = 0.95;
    double lower_quantile = 0.05;

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
    size_t C = dims[1];
    size_t num_samples = dims[0];

    if (condition_a > C || condition_b > C) {
        Logger::abort("Cannot compare conditions %u and %u. The experiment "
                      "has only %u conditions", condition_a, condition_b,
                      C);
    }

    hsize_t mem_dataspace_dims[1] = { T };
    hsize_t mem_dataspace_start[1] = { 0 };
    hid_t mem_dataspace = H5Screate_simple(1, mem_dataspace_dims, NULL);
    H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, mem_dataspace_start,
                        NULL, mem_dataspace_dims, NULL);

    // temporary vector for computing quantiles
    std::vector<float> work(num_samples);

    hsize_t file_dataspace_dims[3] = {1, 1, T};
    hsize_t file_dataspace_start[3] = {0, 0, 0};

    std::vector<float> log2fc_upper(T);
    std::vector<float> log2fc_lower(T);

    std::vector<std::vector<float> >
        condition_data_a(num_samples),
        condition_data_b(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        condition_data_a[i].resize(T);
        condition_data_b[i].resize(T);
    }

    // read ccondition_a data
    for (size_t j = 0; j < num_samples; ++j) {
        file_dataspace_start[0] = j;
        file_dataspace_start[1] = condition_a;

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            file_dataspace_start, NULL,
                            file_dataspace_dims, NULL);

        herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT,
                                mem_dataspace, dataspace, H5P_DEFAULT,
                                &condition_data_a[j].at(0));

        if (status < 0) {
            Logger::abort("Reading /condition/tgroup_mean dataset failed.");
        }
    }

    // read ccondition_b data
    for (size_t j = 0; j < num_samples; ++j) {
        file_dataspace_start[0] = j;
        file_dataspace_start[1] = condition_b;

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            file_dataspace_start, NULL,
                            file_dataspace_dims, NULL);

        herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT,
                                mem_dataspace, dataspace, H5P_DEFAULT,
                                &condition_data_b[j].at(0));

        if (status < 0) {
            Logger::abort("Reading /condition/tgroup_mean dataset failed.");
        }
    }

    for (size_t j = 0; j < T; ++j) {
        for (size_t k = 0; k < num_samples; ++k) {
            work[k] = log2(condition_data_a[k][j]) -
                      log2(condition_data_b[k][j]);
        }

        std::sort(work.begin(), work.end());
        log2fc_lower[j] = work[(int)(lower_quantile * num_samples)];
        log2fc_upper[j] = work[(int)(upper_quantile * num_samples)];
    }

    H5Sclose(mem_dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    fprintf(output, "transcript_ids\tgene_ids\tlower_log2fc\tupper_log2fc\n");
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

        fprintf(output, "\t%f\t%f\n", log2fc_lower[i], log2fc_upper[i]);
    }
}


void Summarize::condition_splicing(std::vector<boost::multi_array<float, 3> >& splicing)
{
    hid_t dataset = H5Dopen2(h5_file, "/condition/splice_mu", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /condition/splice_mu dataset");
    }

    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[3]; // dims are: num_samples, C, spliced_tgroups
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t num_samples = dims[0];
    size_t C = dims[1];

    if (dims[2] != spliced_tgroup_indexes.size()) {
        Logger::abort("The /condition/splice_mu dataset is of incorrect size.");
    }

    hsize_t mem_dataspace_dims[1] = {dims[2]};
    hid_t mem_dataspace = H5Screate_simple(1, mem_dataspace_dims, NULL);

    hvl_t* buffer = new hvl_t[dims[2]];
    hid_t memtype = H5Tvlen_create(H5T_NATIVE_FLOAT);

    hsize_t dataspace_start[3] = {0, 0, 0};
    hsize_t dataspace_count[3] = {1, 1, dims[2]};
    herr_t status;

    for (size_t i = 0; i < dims[2]; ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        splicing[i].resize(boost::extents[dims[0]][dims[1]][tgroup_tids[tgroup].size()]);
    }

    for (size_t i = 0; i < num_samples; ++i) {
        dataspace_start[0] = i;
        for (size_t j = 0; j < C; ++j) {
            dataspace_start[1] = j;
            status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                         dataspace_start, NULL,
                                         dataspace_count, NULL);
            if (status < 0) {
                Logger::abort("HDF5 hyperslab selection failed.");
            }

            status = H5Dread(dataset, memtype, mem_dataspace, dataspace,
                             H5P_DEFAULT, buffer);
            if (status < 0) {
                Logger::abort("HDF5 read failed.");
            }

            for (size_t k = 0; k < dims[2]; ++k) {
                size_t tgroup = spliced_tgroup_indexes[k];
                if (buffer[k].len != tgroup_tids[tgroup].size()) {
                    Logger::abort("Spliced tgroup has an inconsistent transcript count.");
                }

                float sum = 0.0;
                for (size_t l = 0; l < tgroup_tids[tgroup].size(); ++l) {
                    splicing[k][i][j][l] =
                        exp(reinterpret_cast<float*>(buffer[k].p)[l]);
                    sum += splicing[k][i][j][l];
                }

                for (size_t l = 0; l < tgroup_tids[tgroup].size(); ++l) {
                    splicing[k][i][j][l] /= sum;
                }
            }
        }
    }

    H5Dvlen_reclaim(memtype, mem_dataspace, H5P_DEFAULT, buffer);
    delete [] buffer;
    H5Sclose(mem_dataspace);
    H5Tclose(memtype);
    H5Dclose(dataset);
}


void Summarize::condition_splicing(FILE* output)
{
    // TODO: pass these in
    double upper_quantile = 0.95;
    double lower_quantile = 0.05;

    // indexed by sample_num, condition, and within-tgroup transcript
    typedef boost::multi_array<float, 3> marray_t;

    // indexed by spliced tgroup
    std::vector<marray_t> splicing(spliced_tgroup_indexes.size());

    condition_splicing(splicing);
    size_t num_samples = 0;
    size_t C = 0;
    if (splicing.size() > 0) {
        C = splicing[0].shape()[1];
        num_samples = splicing[0].shape()[0];
    }

    // summarization
    fprintf(output, "transcript_id\ttss_transcript_ids");
    for (size_t i = 0; i < C; ++i) {
        unsigned int cond = i + 1;
        fprintf(output, "\tcondition%u_lower\tcondition%u_upper\tcondition%u_median",
                cond, cond, cond);
    }
    fputc('\n', output);

    std::vector<float> work(num_samples);

    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        for (size_t j = 0; j < tgroup_tids[tgroup].size(); ++j) {
            fprintf(output, "%s\t", transcript_ids[tgroup_tids[tgroup][j]].c_str());
            for (size_t k = 0; k < tgroup_tids[tgroup].size(); ++k) {
                if (k != 0) fputc(',', output);
                fprintf(output, "%s", transcript_ids[tgroup_tids[tgroup][k]].c_str());
            }

            for (size_t k = 0; k < C; ++k) {
                for (size_t l = 0; l < num_samples; ++l) {
                    work[l] = splicing[i][l][k][j];
                }

                std::sort(work.begin(), work.end());
                fprintf(output, "\t%f\t%f\t%f",
                        (double) work[lround(work.size() * lower_quantile)],
                        (double) work[lround(work.size() * upper_quantile)],
                        (double) work[work.size()/2]);
            }
            fputc('\n', output);
        }
    }
}


void Summarize::condition_pairwise_splicing(FILE* output)
{
    // TODO: pass these in
    double upper_quantile = 0.95;
    double lower_quantile = 0.05;

    // indexed by sample_num, condition, and within-tgroup transcript
    typedef boost::multi_array<float, 3> marray_t;

    // indexed by spliced tgroup
    std::vector<marray_t> splicing(spliced_tgroup_indexes.size());

    condition_splicing(splicing);
    size_t num_samples = 0;
    size_t C = 0;
    if (splicing.size() > 0) {
        C = splicing[0].shape()[1];
        num_samples = splicing[0].shape()[0];
    }

    fprintf(output, "transcript_id\tgene_id\tcondition_a\tcondition_b\t"
                    "log2fc_lower\tlog2fc_median\tlog2fc_upper\n");

    std::vector<float> work(num_samples);
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        for (size_t j = 0; j < tgroup_tids[tgroup].size(); ++j) {
            unsigned long tid = tgroup_tids[tgroup][j];
            for (size_t cond_a = 0; cond_a < C; ++cond_a) {
                for (size_t cond_b = 0; cond_b < C; ++cond_b) {
                    if (cond_a == cond_b) continue;

                    for (size_t l = 0; l < num_samples; ++l) {
                        work[l] = log2(splicing[i][l][cond_a][j]) -
                                  log2(splicing[i][l][cond_b][j]);
                    }

                    std::sort(work.begin(), work.end());
                    fprintf(output, "%s\t%s\t%lu\t%lu\t%f\t%f\t%f\n",
                            transcript_ids[tid].c_str(), gene_ids[tid].c_str(),
                            cond_a + 1, cond_b + 1,
                            work[(int)(lower_quantile * num_samples)],
                            work[num_samples/2],
                            work[(int)(upper_quantile * num_samples)]);
                }
            }
        }
    }
}


void Summarize::expression_samples(FILE* output)
{
    hid_t dataset = H5Dopen2(h5_file, "/transcript_quantification", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /transcript_quant dataset");
    }

    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[3]; // dims are num_samples, K, N
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t num_samples = dims[0];

    typedef boost::multi_array<float, 3> marray_t;
    marray_t Q(boost::extents[dims[0]][dims[1]][dims[2]]);

    hsize_t file_dataspace_start[3] = {0, 0, 0};
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                        file_dataspace_start, NULL,
                        dims, NULL);
    herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT, dataspace, dataspace,
                            H5P_DEFAULT, &Q.data()[0]);

    if (status < 0) {
        Logger::abort("Reading the /transcript_quantification dataset failed.");
    }

    fprintf(output, "transcript_id");
    for (size_t i = 0; i < K; ++i) {
        fprintf(output, "\tsample%lu", (unsigned long) i);
    }
    fputc('\n', output);

    for (size_t i = 0; i < N; ++i) {
        fputs(transcript_ids[i].c_str(), output);
        for (size_t j = 0; j < K; ++j) {
            fputc('\t', output);
            for (size_t k = 0; k < num_samples; ++k) {
                if (k != 0) fputc(',', output);
                fprintf(output, "%e", (double) Q[k][j][i]);
            }
        }
        fputc('\n', output);
    }

    H5Sclose(dataspace);
    H5Dclose(dataset);
}


void Summarize::read_cassette_exons(std::vector<Interval>& cassette_exons,
                                    std::vector<std::vector<unsigned int> >& including_tids,
                                    std::vector<std::vector<unsigned int> >& excluding_tids)
{
    hid_t dataset = H5Dopen2(h5_file, "/cassette_exons/seqname", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /cassette_exons/seqname dataset");
    }

    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t n = dims[0];

    // read sequence names
    hid_t varstring_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(varstring_type, H5T_VARIABLE);
    std::vector<char*> seqnames(n);
    herr_t status = H5Dread(dataset, varstring_type, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, &seqnames.at(0));
    if (status < 0) {
        Logger::abort("Failed to read the /cassette_exons/seqname dataset.");
    }

    H5Dclose(dataset);

    // read starts
    dataset = H5Dopen2(h5_file, "/cassette_exons/start", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /cassette_exons/start dataset");
    }

    std::vector<long> starts(n);
    status = H5Dread(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &starts.at(0));
    if (status < 0) {
        Logger::abort("Failed to read the /cassette_exons/starts dataset.");
    }

    H5Dclose(dataset);

    // read ends
    dataset = H5Dopen2(h5_file, "/cassette_exons/end", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /cassette_exons/end dataset");
    }

    std::vector<long> ends(n);
    status = H5Dread(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &ends.at(0));
    if (status < 0) {
        Logger::abort("Failed to read the /cassette_exons/ends dataset.");
    }

    H5Dclose(dataset);

    // read strands
    dataset = H5Dopen2(h5_file, "/cassette_exons/strand", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /cassette_exons/strand dataset");
    }

    std::vector<char> strands(n);
    status = H5Dread(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &strands.at(0));
    if (status < 0) {
        Logger::abort("Failed to read the /cassette_exons/strand dataset.");
    }

    H5Dclose(dataset);

    // read including_tids
    dataset = H5Dopen2(h5_file, "/cassette_exons/including_tids", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /cassette_exons/including_tids dataset");
    }

    hid_t tids_type = H5Tvlen_create(H5T_NATIVE_UINT);
    std::vector<hvl_t> including_tids_data(n);
    status = H5Dread(dataset, tids_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &including_tids_data.at(0));
    if (status) {
        Logger::abort("Failed to read the /cassette_exons/including_tids");
    }

    H5Dclose(dataset);

    // read excluding_tids
    dataset = H5Dopen2(h5_file, "/cassette_exons/excluding_tids", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /cassette_exons/excluding_tids dataset");
    }

    std::vector<hvl_t> excluding_tids_data(n);
    status = H5Dread(dataset, tids_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &excluding_tids_data.at(0));
    if (status) {
        Logger::abort("Failed to read the /cassette_exons/excluding_tids");
    }

    H5Dclose(dataset);

    // construct the data
    cassette_exons.resize(n);
    including_tids.resize(n);
    excluding_tids.resize(n);

    for (size_t i = 0; i < n; ++i) {
        cassette_exons[i].seqname = std::string(seqnames[i]);
        cassette_exons[i].start = starts[i];
        cassette_exons[i].end = ends[i];
        if (strands[i] == '+') {
            cassette_exons[i].strand = strand_pos;
        }
        else if (strands[i] == '-') {
            cassette_exons[i].strand = strand_neg;
        }
        else {
            cassette_exons[i].strand = strand_na;
        }

        including_tids[i].resize(including_tids_data[i].len);
        for (size_t j = 0; j < including_tids_data[i].len; ++j) {
            including_tids[i][j] =
                reinterpret_cast<unsigned int*>(including_tids_data[i].p)[j];
        }

        excluding_tids[i].resize(excluding_tids_data[i].len);
        for (size_t j = 0; j < excluding_tids_data[i].len; ++j) {
            excluding_tids[i][j] =
                reinterpret_cast<unsigned int*>(excluding_tids_data[i].p)[j];
        }
    }

    H5Dvlen_reclaim(varstring_type, dataspace, H5P_DEFAULT, &seqnames.at(0));
    H5Tclose(varstring_type);

    H5Dvlen_reclaim(tids_type, dataspace, H5P_DEFAULT, &including_tids_data.at(0));
    H5Dvlen_reclaim(tids_type, dataspace, H5P_DEFAULT, &excluding_tids_data.at(0));
    H5Tclose(tids_type);
}


void Summarize::read_tgroup_mean(boost::multi_array<float, 3>& output)
{
    hid_t dataset = H5Dopen2(h5_file, "/condition/tgroup_mean", H5P_DEFAULT);
    if (dataset < 0) {
        Logger::abort("Failed to open the /condition/tgroup_mean dataset");
    }

    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[3]; // dims are: num_samples, C, T
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    output.resize(boost::extents[dims[0]][dims[1]][dims[2]]);

    herr_t status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, output.data());
    if (status < 0) {
        Logger::abort("HDF5 read failed.");
    }

    H5Dclose(dataset);
}


void Summarize::cassette_exon_pairwise_splicing(FILE* output)
{
    // TODO: pass these in
    double upper_quantile = 0.95;
    double lower_quantile = 0.05;

    // indexed by sample_num, condition, and within-tgroup transcript
    typedef boost::multi_array<float, 3> marray_t;

    // indexed by spliced tgroup
    std::vector<marray_t> splicing(spliced_tgroup_indexes.size());

    condition_splicing(splicing);
    size_t num_samples = 0;
    size_t C = 0;
    if (splicing.size() > 0) {
        C = splicing[0].shape()[1];
        num_samples = splicing[0].shape()[0];
    }

    boost::multi_array<float, 3> tgroup_means;
    read_tgroup_mean(tgroup_means);

    std::vector<Interval> cassette_exons;
    std::vector<std::vector<unsigned int> > including_tids;
    std::vector<std::vector<unsigned int> > excluding_tids;

    read_cassette_exons(cassette_exons, including_tids, excluding_tids);


    // construct an index mapping tids to their index within their tgroup
    std::vector<unsigned int> tid_tgroup_index(N);
    BOOST_FOREACH (const std::vector<unsigned int>& tids, tgroup_tids) {
        for (size_t i = 0; i < tids.size(); ++i) {
            tid_tgroup_index[tids[i]] = i;
        }
    }

    // construct an index mapping tgroup to a spliced tgroup index
    std::vector<unsigned int> tgroup_spliced_tgroup(T);
    std::fill(tgroup_spliced_tgroup.begin(), tgroup_spliced_tgroup.end(), (unsigned int) -1);
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        tgroup_spliced_tgroup[spliced_tgroup_indexes[i]] = i;
    }

    fprintf(output,
            "seqname\tstart\tend\tstrand\tspliced_in_transcript_ids\t"
            "spliced_out_transcript_ids\t"
            "condition_a\tcondition_b\tlog2fc_lower\t"
            "log2fc_median\tlog2fc_upper\n");

    boost::multi_array<float, 2> condition_spliced_in_proportion(
            boost::extents[C][num_samples]);
    std::vector<float> work(num_samples);

    for (size_t i = 0; i < cassette_exons.size(); ++i) {
        for (size_t j = 0; j < C; ++j) {
            for (size_t k = 0; k < num_samples; ++k) {
                float spliced_in = 0.0;
                BOOST_FOREACH (unsigned int tid, including_tids[i]) {
                    unsigned int tg = tgroup[tid];
                    float tmix = 0.0;
                    if (tgroup_spliced_tgroup[tg] < (unsigned int) -1) {
                        tmix = log(splicing[tgroup_spliced_tgroup[tg]][k][j][tid_tgroup_index[tid]]);
                    }

                    spliced_in += tmix + tgroup_means[k][j][tg];
                }

                float spliced_out = 0.0;
                BOOST_FOREACH (unsigned int tid, excluding_tids[i]) {
                    unsigned int tg = tgroup[tid];
                    float tmix = 0.0;
                    if (tgroup_spliced_tgroup[tg] < (unsigned int) -1) {
                        tmix = log(splicing[tgroup_spliced_tgroup[tg]][k][j][tid_tgroup_index[tid]]);
                    }

                    spliced_out += tmix + tgroup_means[k][j][tg];
                }

                // compute the proportion spliced in
                condition_spliced_in_proportion[j][k] =
                    (spliced_in - logaddexp(spliced_in, spliced_out)) / M_LN2;
            }
        }

        for (size_t cond_a = 0; cond_a < C; ++cond_a) {
            for (size_t cond_b = 0; cond_b < C; ++cond_b) {
                if (cond_a == cond_b) continue;

                for (size_t k = 0; k < num_samples; ++k) {
                    work[k] = condition_spliced_in_proportion[cond_a][k] -
                              condition_spliced_in_proportion[cond_b][k];
                }

                std::sort(work.begin(), work.end());

                fprintf(output, "%s\t%ld\t%ld\t",
                        cassette_exons[i].seqname.get().c_str(),
                        cassette_exons[i].start,
                        cassette_exons[i].end);

                switch (cassette_exons[i].strand) {
                    case strand_pos:
                        fputc('+', output);
                        break;
                    case strand_neg:
                        fputc('-', output);
                        break;
                    default:
                        fputc('.', output);
                }
                fputc('\t', output);


                bool first = true;
                BOOST_FOREACH (unsigned int tid, including_tids[i]) {
                    if (!first) {
                        fputc(',', output);
                    }
                    else first = false;
                    fputs(transcript_ids[tid].c_str(), output);
                }
                fputc('\t', output);

                first = true;
                BOOST_FOREACH (unsigned int tid, excluding_tids[i]) {
                    if (!first) {
                        fputc(',', output);
                    }
                    else first = false;
                    fputs(transcript_ids[tid].c_str(), output);
                }
                fputc('\t', output);

                fprintf(output, "%lu\t%lu\t%f\t%f\t%f\n",
                        cond_a + 1, cond_b + 1,
                        work[(int)(lower_quantile * num_samples)],
                        work[num_samples/2],
                        work[(int)(upper_quantile * num_samples)]);
            }
        }
    }
}


Summarize::~Summarize()
{
    H5Fclose(h5_file);
}

