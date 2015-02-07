
#include <algorithm>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <intervals.hpp>
#include <set>

#include "constants.hpp"
#include "summarize.hpp"
#include "logger.hpp"


using namespace boost::numeric::ublas;


// Compute log(exp(x) + exp(y)), avoiding overflow/underflow.
static double logaddexp(double x, double y)
{
    double u = x - y;
    if (u > 0.0) {
        return x + log1p(exp(-u));
    } else if (u <= 0.0) {
        return y + log1p(exp(u));
    } else  {
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

    dataset = H5Dopen2_checked(h5_file, "/transcript_id", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    N = dims[0];

    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);

    const char** string_data = new const char* [N];

    // read transcript_id dataset
    H5Dread_checked(dataset, datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, string_data);

    transcript_ids.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        transcript_ids.push_back(TranscriptID(std::string(string_data[i])));
    }

    H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, string_data);
    H5Dclose(dataset);

    // read gene_name dataset
    dataset = H5Dopen2_checked(h5_file, "/gene_name", H5P_DEFAULT);
    H5Dread_checked(dataset, datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, string_data);

    for (size_t i = 0; i < N; ++i) {
        gene_names.push_back(GeneName(std::string(string_data[i])));
    }

    H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, string_data);
    H5Dclose(dataset);

    // read gene_id dataset
    dataset = H5Dopen2_checked(h5_file, "/gene_id", H5P_DEFAULT);
    H5Dread_checked(dataset, datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, string_data);

    gene_ids.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        gene_ids.push_back(GeneID(std::string(string_data[i])));
    }

    H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT, string_data);
    H5Dclose(dataset);

    delete [] string_data;
    H5Tclose(datatype);

    for (size_t i = 0; i < N; ++i) {
        gid_to_tids[gene_ids[i]].push_back(i);
        gid_to_gene_name[gene_ids[i]] = gene_names[i];
    }

    // read tgroups
    dataset = H5Dopen2_checked(h5_file, "/tgroup", H5P_DEFAULT);

    tgroup.resize(N);
    H5Dread_checked(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, &tgroup.at(0));

    T = *std::max_element(tgroup.begin(), tgroup.end()) + 1;

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
    dataset = H5Dopen2_checked(h5_file, "/transcript_quantification", H5P_DEFAULT);

    dataspace = H5Dget_space(dataset);
    hsize_t dims3[3]; // dims are num_samples, K, N
    H5Sget_simple_extent_dims(dataspace, dims3, NULL);
    K = dims3[1];

    H5Dclose(dataset);
    H5Sclose(dataspace);

    num_samples = dims3[0];

    // read scale
    dataset = H5Dopen2_checked(h5_file, "/sample_scaling", H5P_DEFAULT);

    scale.resize(num_samples, K);
    H5Dread_checked(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, &scale.data()[0]);

    H5Dclose(dataset);

    read_metadata(metadata);

    std::string last_condition_name;
    BOOST_FOREACH (std::string& condition_name, metadata.sample_conditions) {
        if (condition_name != last_condition_name) {
            condition_names.push_back(condition_name);
            last_condition_name = condition_name;
        }
    }

    C = condition_names.size();
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

    hid_t dataset = H5Dopen2_checked(h5_file, "/condition/tgroup_mean", H5P_DEFAULT);

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

            H5Dread_checked(dataset, H5T_NATIVE_FLOAT,
                            mem_dataspace, dataspace, H5P_DEFAULT,
                            &condition_data[j].at(0));
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

    hid_t dataset = H5Dopen2_checked(h5_file, "/experiment/tgroup_sd", H5P_DEFAULT);

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

        H5Dread_checked(dataset, H5T_NATIVE_FLOAT, mem_dataspace,
                        dataspace, H5P_DEFAULT, &data[i].at(0));
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


void Summarize::point_ci_transcript_expression(
        matrix<float>* point, matrix<float>* lower, matrix<float>* upper,
        double interval, bool unnormalized, bool splicing_rate)
{
    double lower_quantile = 0.5 - interval/2;
    double upper_quantile = 0.5 + interval/2;

    hid_t dataset = H5Dopen2_checked(h5_file, "/transcript_quantification", H5P_DEFAULT);

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
        file_dataspace_start[1] = i;

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            file_dataspace_start, NULL,
                            file_dataspace_dims, NULL);

        H5Dread_checked(dataset, H5T_NATIVE_FLOAT, mem_dataspace,
                        dataspace, H5P_DEFAULT, &Qi.data()[0]);

        if (unnormalized) {
            for (size_t j = 0; j < num_samples; ++j) {
                for (size_t k = 0; k < N; ++k) {
                    Qi(j, k) /= scale(j, i);
                }
            }
        }

        // truncate miniscule values
        for (size_t j = 0; j < num_samples; ++j) {
            for (size_t k = 0; k < N; ++k) {
                Qi(j, k) = std::max<float>(constants::min_expr, Qi(j, k));
            }
        }

        if (splicing_rate) {
            for (size_t tg = 0; tg < T; ++tg) {
                for (size_t j = 0; j < num_samples; ++j) {
                    double expr_sum = 0.0;
                    BOOST_FOREACH (unsigned int tid, tgroup_tids[tg]) {
                        expr_sum += Qi(j, tid);
                    }

                    BOOST_FOREACH (unsigned int tid, tgroup_tids[tg]) {
                        Qi(j, tid) /= expr_sum;
                    }

                    // XXX: experiment with seeing the transformed data
                    //BOOST_FOREACH (unsigned int tid, tgroup_tids[tg]) {
                        //Qi(j, tid) =
                            //log(Qi(j, tid) / Qi(j, tgroup_tids[tg].back()));
                    //}
                }
            }
        }

        for (size_t j = 0; j < N; ++j) {
            matrix_column<matrix<float> > col(Qi, j);

            // TODO: maximum posterior. there should be a switch to use this
            //(*point)(i, j) = col[0];

            // posterior mean
            std::sort(col.begin(), col.end());
            (*point)(i, j) = col[col.size() / 2];

            if (lower) {
                (*lower)(i, j) = col[lround((col.size() - 1) * lower_quantile)];
                if ((*point)(i, j) < (*lower)(i, j)) (*lower)(i, j) = (*point)(i, j);
            }
            if (upper) {
                (*upper)(i, j) = col[lround((col.size() - 1) * upper_quantile)];
                if ((*point)(i, j) > (*upper)(i, j)) (*upper)(i, j) = (*point)(i, j);
            }
        }
    }

    H5Sclose(mem_dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}


void Summarize::point_ci_gene_expression(
        matrix<float>* point, matrix<float>* lower, matrix<float>* upper,
        double interval, bool unnormalized)
{
    double lower_quantile = 0.5 - interval/2;
    double upper_quantile = 0.5 + interval/2;

    hid_t dataset = H5Dopen2_checked(h5_file, "/transcript_quantification", H5P_DEFAULT);

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
    std::vector<float> work(num_samples);

    typedef std::pair<std::string, std::vector<size_t> > item_t;
    for (size_t i = 0; i < K; ++i) {
        file_dataspace_start[1] = i;

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            file_dataspace_start, NULL,
                            file_dataspace_dims, NULL);

        H5Dread_checked(dataset, H5T_NATIVE_FLOAT, mem_dataspace,
                        dataspace, H5P_DEFAULT, &Qi.data()[0]);

        if (unnormalized) {
            for (size_t j = 0; j < num_samples; ++j) {
                for (size_t k = 0; k < N; ++k) {
                    Qi(j, k) /= scale(j, i);
                }
            }
        }

        size_t j = 0;
        BOOST_FOREACH (const item_t& item, gid_to_tids) {
            std::fill(work.begin(), work.end(), 0.0);
            BOOST_FOREACH (size_t tid, item.second) {
                for (size_t k = 0; k < num_samples; ++k) {
                    work[k] += Qi(k, tid);
                }
            }

            // truncate miniscule values
            for (size_t k = 0; k < num_samples; ++k) {
                work[k] = std::max<float>(constants::min_expr, work[k]);
            }

            std::sort(work.begin(), work.end());

            (*point)(i, j) = work[work.size() / 2];

            if (lower) {
                (*lower)(i, j) = work[lround((num_samples - 1) * lower_quantile)];
                if ((*point)(i, j) < (*lower)(i, j)) (*lower)(i, j) = (*point)(i, j);
            }
            if (upper) {
                (*upper)(i, j) = work[lround((num_samples - 1) * upper_quantile)];
                if ((*point)(i, j) > (*upper)(i, j)) (*upper)(i, j) = (*point)(i, j);
            }
            ++j;
        }
    }

    H5Sclose(mem_dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}


void Summarize::transcript_expression(FILE* output, double credible_interval,
                                      bool unnormalized, bool splicing_rate)
{
    bool print_credible_interval = !isnan(credible_interval);

    matrix<float> point(K, N);
    matrix<float> lower;
    matrix<float> upper;

    if (print_credible_interval) {
        lower.resize(K, N);
        upper.resize(K, N);
        point_ci_transcript_expression(&point, &lower, &upper,
                                        credible_interval, unnormalized,
                                        splicing_rate);
    }
    else {
        point_ci_transcript_expression(&point, NULL, NULL, credible_interval,
                                        unnormalized, splicing_rate);
    }

    fprintf(output, "gene_name\tgene_id\ttranscript_id");
    for (unsigned int i = 0; i < K; ++i) {
        fprintf(output,
                splicing_rate ? "\t%s_splice_rate" :
                    (unnormalized ?  "\t%s_tpm" : "\t%s_adjusted_tpm"),
                metadata.sample_names[i].c_str());
        if (print_credible_interval) {
            fprintf(output,
                    splicing_rate ? "\t%s_splice_rate_lower\t%s_splice_rate_upper" :
                    (unnormalized ?
                        "\t%s_tpm_lower\t%s_tpm_upper" :
                        "\t%s_adjusted_tpm_lower\t%s_adjusted_tpm_upper"),
                    metadata.sample_names[i].c_str(),
                    metadata.sample_names[i].c_str());
        }
    }
    fputc('\n', output);

    double expr_scale = splicing_rate ? 1.0 : 1e6;

    for (size_t i = 0; i < N; ++i) {
        if (splicing_rate && tgroup_tids[tgroup[i]].size() < 2) {
            continue;
        }

        fprintf(output, "%s\t%s\t%s",
                gene_names[i].get().c_str(),
                gene_ids[i].get().c_str(),
                transcript_ids[i].get().c_str());

        for (size_t j = 0; j < K; ++j) {
            if (print_credible_interval) {
                fprintf(output, "\t%e\t%e\t%e",
                        expr_scale * point(j, i),
                        expr_scale * lower(j, i),
                        expr_scale * upper(j, i));
            }
            else {
                fprintf(output, "\t%e", expr_scale * point(j, i));
            }
        }
        fputc('\n', output);
    }
}


void Summarize::gene_expression(FILE* output, double credible_interval,
                                bool unnormalized)
{
    bool print_credible_interval = !isnan(credible_interval);
    size_t num_genes = gid_to_tids.size();

    matrix<float> point(K, num_genes);
    matrix<float> lower;
    matrix<float> upper;

    if (print_credible_interval) {
        lower.resize(K, num_genes);
        upper.resize(K, num_genes);
        point_ci_gene_expression(&point, &lower, &upper,
                                  credible_interval, unnormalized);
    }
    else {
        point_ci_gene_expression(&point, NULL, NULL, credible_interval,
                                  unnormalized);
    }

    fprintf(output, "gene_name\tgene_id\ttranscript_ids");
    for (unsigned int i = 0; i < K; ++i) {
        fprintf(output,
                unnormalized ? "\t%s_tpm" : "\t%s_adjusted_tpm",
                metadata.sample_names[i].c_str());
        if (print_credible_interval) {
            fprintf(output,
                    unnormalized ?
                        "\t%s_tpm_lower\t%s_tpm_upper" :
                        "\t%s_adjusted_tpm_lower\t%s_adjusted_tpm_upper",
                    metadata.sample_names[i].c_str(),
                    metadata.sample_names[i].c_str());
        }
    }
    fputc('\n', output);

    size_t i = 0;
    typedef std::pair<GeneID, std::vector<size_t> > item_t;
    BOOST_FOREACH (const item_t& gid_tids, gid_to_tids) {
        const GeneID& gene_id = gid_tids.first;
        const std::vector<size_t>& tids = gid_tids.second;
        const GeneName& gene_name = gid_to_gene_name[gene_id];

        fprintf(output, "%s\t%s\t%s",
                gene_name.get().c_str(),
                gene_id.get().c_str(),
                transcript_ids[tids[0]].get().c_str());

        for (size_t j = 1; j < tids.size(); ++j) {
            fprintf(output, ",%s", transcript_ids[tids[j]].get().c_str());
        }

        for (size_t j = 0; j < K; ++j) {
            if (print_credible_interval) {
                fprintf(output, "\t%e\t%e\t%e",
                        1e6 * point(j, i), 1e6 * lower(j, i), 1e6 * upper(j, i));
            }
            else {
                fprintf(output, "\t%e", 1e6 * point(j, i));
            }
        }
        fputc('\n', output);
        ++i;
    }
}


void Summarize::read_condition_tgroup_mean(unsigned int condition,
                                           matrix<float>& data)
{
    hid_t dataset = H5Dopen2_checked(h5_file, "/condition/tgroup_mean",
                                     H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    hsize_t mem_dataspace_start[2] = {0, 0};
    hsize_t mem_dataspace_dims[2] = {num_samples, T};
    hid_t mem_dataspace = H5Screate_simple(2, mem_dataspace_dims, NULL);
    H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, mem_dataspace_start,
                        NULL, mem_dataspace_dims, NULL);

    hsize_t file_dataspace_start[3] = {0, condition, 0};
    hsize_t file_dataspace_dims[3] = {num_samples, 1, T};
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, file_dataspace_start, NULL,
                        file_dataspace_dims, NULL);
    H5Dread_checked(dataset, H5T_NATIVE_FLOAT, mem_dataspace, dataspace,
                    H5P_DEFAULT, &data.data()[0]);

    H5Sclose(mem_dataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}



void Summarize::differential_transcription(FILE* output, double credible_interval,
                                           double effect_size)
{
    if (isnan(effect_size)) {
        effect_size = 1.0;
    }
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    fprintf(output,
            "gene_names\tgene_ids\ttranscript_ids\tcondition_a\tcondition_b\t"
            "down_pr\tup_pr\tmedian_log2_fold_change");
    if (print_credible_interval) {
        fprintf(output, "\tlower_log2_fold_change\tupper_log2_fold_change");
    }
    fputc('\n', output);

    hid_t dataset = H5Dopen2_checked(h5_file, "/condition/tgroup_mean",
                                     H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[3]; // dims are: num_samples, C, spliced_tgroups
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    matrix<float> condition_a_data(num_samples, T);
    matrix<float> condition_b_data(num_samples, T);
    std::vector<float> work(num_samples);

    size_t C = dims[1];
    for (unsigned int condition_a = 0; condition_a < C - 1; ++condition_a) {
        read_condition_tgroup_mean(condition_a, condition_a_data);
        for (unsigned int condition_b = condition_a + 1; condition_b < C; ++condition_b) {
            read_condition_tgroup_mean(condition_b, condition_b_data);

            for (size_t i = 0; i < T; ++i) {
                // gene names
                std::set<GeneName> tgroup_gene_names;
                BOOST_FOREACH (unsigned int tid, tgroup_tids[i]) {
                    tgroup_gene_names.insert(gene_names[tid]);
                }
                bool first = true;
                BOOST_FOREACH (const GeneName& gene_name, tgroup_gene_names) {
                    if (first) first = false;
                    else fputc(',', output);
                    fputs(gene_name.get().c_str(), output);
                }

                // gene ids
                fputc('\t', output);
                std::set<GeneID> tgroup_gene_ids;
                BOOST_FOREACH (unsigned int tid, tgroup_tids[i]) {
                    tgroup_gene_ids.insert(gene_ids[tid]);
                }
                first = true;
                BOOST_FOREACH (const GeneID& gene_id, tgroup_gene_ids) {
                    if (first) first = false;
                    else fputc(',', output);
                    fputs(gene_id.get().c_str(), output);
                }

                // transcript ids
                fputc('\t', output);
                first = true;
                BOOST_FOREACH (unsigned int tid, tgroup_tids[i]) {
                    if (first) first = false;
                    else fputc(',', output);
                    fputs(transcript_ids[tid].get().c_str(), output);
                }

                // conditions
                fprintf(output, "\t%s\t%s",
                        condition_names[condition_a].c_str(),
                        condition_names[condition_b].c_str());

                double down_pr = 0, up_pr = 0;
                for (size_t j = 0; j < num_samples; ++j) {
                    work[j] = (condition_a_data(j, i) -
                               condition_b_data(j, i)) / M_LN2;
                }
                std::sort(work.begin(), work.end());

                BOOST_FOREACH (float l2fc, work) {
                    if (fabs(l2fc) >= effect_size) {
                        if (l2fc <= 0.0) down_pr += 1;
                        else             up_pr += 1;
                    }
                }
                down_pr /= num_samples;
                up_pr /= num_samples;
                fprintf(output, "\t%0.3f\t%0.3f\t%e",
                        down_pr, up_pr, work[num_samples/2]);

                if (print_credible_interval) {
                    fprintf(output, "\t%e\t%e",
                            work[lround((num_samples - 1) * lower_quantile)],
                            work[lround((num_samples - 1) * upper_quantile)]);
                }
                fputc('\n', output);
            }
        }
    }
}


void Summarize::differential_splicing(FILE* output, double credible_interval,
                                      double effect_size)
{
    if (isnan(effect_size)) {
        effect_size = 0.2;
    }
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

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

    std::vector<float> work(num_samples);

    fprintf(output,
            "gene_name\tgene_id\ttranscript_id\tcondition_a\tcondition_b\t"
            "down_pr\tup_pr\tmedian_change");
    if (print_credible_interval) {
        fprintf(output, "\tlower_change\tupper_change");
    }
    fputc('\n', output);

    for (unsigned int condition_a = 0; condition_a < C - 1; ++condition_a) {
        for (unsigned int condition_b = condition_a + 1; condition_b < C; ++condition_b) {
            for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
                size_t tgroup = spliced_tgroup_indexes[i];
                for (size_t j = 0; j < tgroup_tids[tgroup].size(); ++j) {
                    unsigned long tid = tgroup_tids[tgroup][j];
                    double down_pr = 0, up_pr = 0;
                    for (size_t l = 0; l < num_samples; ++l) {
                        work[l] = splicing[i][l][condition_a][j] -
                                  splicing[i][l][condition_b][j];
                        if (fabs(work[l]) >= effect_size) {
                            if (work[l] <= 0) down_pr += 1;
                            else up_pr += 1;
                        }
                    }
                    down_pr /= num_samples;
                    up_pr /= num_samples;

                    std::sort(work.begin(), work.end());
                    fprintf(output, "%s\t%s\t%s\t%s\t%s\t%0.3f\t%0.3f\t%e",
                            gene_names[tid].get().c_str(),
                            gene_ids[tid].get().c_str(),
                            transcript_ids[tid].get().c_str(),
                            condition_names[condition_a].c_str(),
                            condition_names[condition_b].c_str(),
                            down_pr, up_pr, work[num_samples/2]);
                    if (print_credible_interval) {
                        fprintf(output, "\t%e\t%e",
                                work[lround((num_samples - 1) * lower_quantile)],
                                work[lround((num_samples - 1) * upper_quantile)]);
                    }
                    fputc('\n', output);
                }
            }
        }
    }
}


void Summarize::condition_splicing(std::vector<boost::multi_array<float, 3> >& splicing)
{
    hid_t dataset = H5Dopen2_checked(h5_file, "/condition/splice_mu", H5P_DEFAULT);

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

            H5Dread_checked(dataset, memtype, mem_dataspace, dataspace,
                            H5P_DEFAULT, buffer);

            for (size_t k = 0; k < dims[2]; ++k) {
                size_t tgroup = spliced_tgroup_indexes[k];
                if (buffer[k].len != tgroup_tids[tgroup].size()) {
                    Logger::abort("Spliced tgroup has an inconsistent transcript count.");
                }

                for (size_t l = 0; l < tgroup_tids[tgroup].size(); ++l) {
                    splicing[k][i][j][l] =
                        std::min<double>(1.0,
                            std::max<double>(0.0,
                                reinterpret_cast<float*>(buffer[k].p)[l]));
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


void Summarize::experiment_splicing(FILE* output, double credible_interval)
{
    // TODO: most of this code is identical to experiment_splicing_sigma

    hid_t dataset = H5Dopen2_checked(h5_file, "/experiment/splice_mu", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    hsize_t dims[3]; // dims are: num_samples, spliced_tgroups
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t num_samples = dims[0];

    if (dims[1] != spliced_tgroup_indexes.size()) {
        Logger::abort("The /experiment/splice_sigma dataset is of incorrect size.");
    }

    hsize_t mem_dataspace_dims[1] = {dims[1]};
    hid_t mem_dataspace = H5Screate_simple(1, mem_dataspace_dims, NULL);

    hvl_t* buffer = new hvl_t[dims[1]];
    hid_t memtype = H5Tvlen_create(H5T_NATIVE_FLOAT);

    hsize_t dataspace_start[3] = {0, 0};
    hsize_t dataspace_count[3] = {1, dims[1]};

    std::vector<boost::multi_array<float, 2> > mu(spliced_tgroup_indexes.size());
    for (size_t i = 0; i < dims[1]; ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        mu[i].resize(boost::extents[dims[0]][tgroup_tids[tgroup].size()]);
    }

    for (size_t i = 0; i < num_samples; ++i) {
        dataspace_start[0] = i;

        H5Sselect_hyperslab_checked(dataspace, H5S_SELECT_SET,
                                    dataspace_start, NULL,
                                    dataspace_count, NULL);

        H5Dread_checked(dataset, memtype, mem_dataspace, dataspace,
                        H5P_DEFAULT, buffer);

        for (size_t k = 0; k < dims[1]; ++k) {
            size_t tgroup = spliced_tgroup_indexes[k];
            if (buffer[k].len != tgroup_tids[tgroup].size()) {
                Logger::abort("Spliced tgroup has an inconsistent transcript count.");
            }

            for (size_t l = 0; l < tgroup_tids[tgroup].size(); ++l) {
                mu[k][i][l] = reinterpret_cast<float*>(buffer[k].p)[l];
            }

            for (size_t l = 0; l < tgroup_tids[tgroup].size(); ++l) {
                mu[k][i][l] = reinterpret_cast<float*>(buffer[k].p)[l];
            }
        }
    }
    H5Dvlen_reclaim(memtype, mem_dataspace, H5P_DEFAULT, buffer);
    delete [] buffer;
    H5Sclose(mem_dataspace);
    H5Tclose(memtype);
    H5Dclose(dataset);

    // print stuff
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    fprintf(output, "gene_names\tgene_ids\ttranscription_group\ttranscript_id");
    if (print_credible_interval) {
        fprintf(output, "\tsplice");
    }
    else {
        fprintf(output, "\tsplice\tsplice_lower\tsplice_upper");
    }
    fputc('\n', output);

    std::set<GeneID> tgroup_gene_ids;
    std::set<GeneName> tgroup_gene_names;

    std::vector<float> work(num_samples);
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        tgroup_gene_names.clear();
        tgroup_gene_ids.clear();
        BOOST_FOREACH (unsigned int tid, tgroup_tids[tgroup]) {
            tgroup_gene_names.insert(gene_names[tid]);
            tgroup_gene_ids.insert(gene_ids[tid]);
        }

        for (size_t j = 0; j < tgroup_tids[tgroup].size(); ++j) {
            unsigned int tid = tgroup_tids[tgroup][j];

            // gene_names
            bool first_item = true;
            BOOST_FOREACH (const GeneName& gene_name, tgroup_gene_names) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", gene_name.get().c_str());
            }
            fputc('\t', output);

            // gene_ids
            first_item = true;
            BOOST_FOREACH (const GeneID& gene_id, tgroup_gene_ids) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", gene_id.get().c_str());
            }
            fputc('\t', output);

            // transcription_group
            first_item = true;
            BOOST_FOREACH (unsigned int tid, tgroup_tids[tgroup]) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", transcript_ids[tid].get().c_str());
            }
            fputc('\t', output);

            // transcript_id
            fprintf(output, "%s", transcript_ids[tid].get().c_str());

            for (size_t l = 0; l < num_samples; ++l) {
                work[l] = mu[i][l][j];
            }
            std::sort(work.begin(), work.end());

            if (print_credible_interval) {
                fprintf(output, "\t%f\t%f\t%f",
                        (double) work[lround((num_samples -1) * lower_quantile)],
                        (double) work[num_samples/2],
                        (double) work[lround((num_samples -1) * upper_quantile)]);
            }
            else {
                fprintf(output, "\t%f", (double) work[num_samples/2]);
            }
            fputc('\n', output);
        }
    }
}


void Summarize::experiment_splicing_sigma(FILE* output, double credible_interval)
{
    hid_t dataset = H5Dopen2_checked(h5_file, "/experiment/splice_sigma", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    hsize_t dims[3]; // dims are: num_samples, spliced_tgroups
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t num_samples = dims[0];

    if (dims[1] != spliced_tgroup_indexes.size()) {
        Logger::abort("The /experiment/splice_sigma dataset is of incorrect size.");
    }

    hsize_t mem_dataspace_dims[1] = {dims[1]};
    hid_t mem_dataspace = H5Screate_simple(1, mem_dataspace_dims, NULL);

    hvl_t* buffer = new hvl_t[dims[1]];
    hid_t memtype = H5Tvlen_create(H5T_NATIVE_FLOAT);

    hsize_t dataspace_start[3] = {0, 0};
    hsize_t dataspace_count[3] = {1, dims[1]};

    std::vector<boost::multi_array<float, 2> > sigma(spliced_tgroup_indexes.size());
    for (size_t i = 0; i < dims[1]; ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        sigma[i].resize(boost::extents[dims[0]][tgroup_tids[tgroup].size()]);
    }

    for (size_t i = 0; i < num_samples; ++i) {
        dataspace_start[0] = i;

        H5Sselect_hyperslab_checked(dataspace, H5S_SELECT_SET,
                                    dataspace_start, NULL,
                                    dataspace_count, NULL);

        H5Dread_checked(dataset, memtype, mem_dataspace, dataspace,
                        H5P_DEFAULT, buffer);

        for (size_t k = 0; k < dims[1]; ++k) {
            size_t tgroup = spliced_tgroup_indexes[k];
            if (buffer[k].len != tgroup_tids[tgroup].size()) {
                Logger::abort("Spliced tgroup has an inconsistent transcript count.");
            }
            for (size_t l = 0; l < tgroup_tids[tgroup].size(); ++l) {
                sigma[k][i][l] = reinterpret_cast<float*>(buffer[k].p)[l];
            }
        }
    }
    H5Dvlen_reclaim(memtype, mem_dataspace, H5P_DEFAULT, buffer);
    delete [] buffer;
    H5Sclose(mem_dataspace);
    H5Tclose(memtype);
    H5Dclose(dataset);

    // print stuff
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    fprintf(output, "gene_names\tgene_ids\ttranscription_group\ttranscript_id");
    if (print_credible_interval) {
        fprintf(output, "\tsigma");
    }
    else {
        fprintf(output, "\tsigma\tsigma_lower\tsigma_upper");
    }
    fputc('\n', output);

    std::set<GeneID> tgroup_gene_ids;
    std::set<GeneName> tgroup_gene_names;

    std::vector<float> work(num_samples);
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        tgroup_gene_names.clear();
        tgroup_gene_ids.clear();
        BOOST_FOREACH (unsigned int tid, tgroup_tids[tgroup]) {
            tgroup_gene_names.insert(gene_names[tid]);
            tgroup_gene_ids.insert(gene_ids[tid]);
        }

        for (size_t j = 0; j < tgroup_tids[tgroup].size(); ++j) {
            unsigned int tid = tgroup_tids[tgroup][j];

            // gene_names
            bool first_item = true;
            BOOST_FOREACH (const GeneName& gene_name, tgroup_gene_names) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", gene_name.get().c_str());
            }
            fputc('\t', output);

            // gene_ids
            first_item = true;
            BOOST_FOREACH (const GeneID& gene_id, tgroup_gene_ids) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", gene_id.get().c_str());
            }
            fputc('\t', output);

            // transcription_group
            first_item = true;
            BOOST_FOREACH (unsigned int tid, tgroup_tids[tgroup]) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", transcript_ids[tid].get().c_str());
            }
            fputc('\t', output);

            // transcript_id
            fprintf(output, "%s", transcript_ids[tid].get().c_str());

            for (size_t l = 0; l < num_samples; ++l) {
                work[l] = sigma[i][l][j];
            }
            std::sort(work.begin(), work.end());

            if (print_credible_interval) {
                fprintf(output, "\t%f\t%f\t%f",
                        (double) work[lround((num_samples -1) * lower_quantile)],
                        (double) work[num_samples/2],
                        (double) work[lround((num_samples -1) * upper_quantile)]);
            }
            else {
                fprintf(output, "\t%f", (double) work[num_samples/2]);
            }
            fputc('\n', output);
        }
    }
}


void Summarize::condition_splicing_sigma(FILE* output, double credible_interval)
{
    // TODO: implement this function
    UNUSED(output);
    UNUSED(credible_interval);
}


std::vector<float> Summarize::transcript_experiment_splicing(const char* transcript_id)
{
    std::vector<float> output(num_samples);

    size_t tid = 0;
    for (; tid < transcript_ids.size(); ++tid) {
        if (transcript_ids[tid] == transcript_id) break;
    }

    if (tid == transcript_ids.size()) {
        Logger::abort("Transcript ID %s not present in dataset.", transcript_id);
    }

    size_t tg = tgroup[tid];

    size_t stg = 0;
    for (; stg < spliced_tgroup_indexes.size(); ++stg) {
        if (spliced_tgroup_indexes[stg] == tg) break;
    }

    if (stg == spliced_tgroup_indexes.size()) {
        Logger::abort("%s is not alternatively spliced.", transcript_id);
    }

    size_t idx = 0;
    for (; idx < tgroup_tids[tg].size(); ++idx) {
        if (tgroup_tids[tg][idx] == tid) break;
    }

    hid_t dataset = H5Dopen2_checked(h5_file, "/experiment/splice_mu", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    hsize_t dims[3]; // dims are: num_samples, spliced_tgroups
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t num_samples = dims[0];

    if (dims[1] != spliced_tgroup_indexes.size()) {
        Logger::abort("The /experiment/splice_sigma dataset is of incorrect size.");
    }

    hsize_t mem_dataspace_dims[1] = {dims[1]};
    hid_t mem_dataspace = H5Screate_simple(1, mem_dataspace_dims, NULL);

    hvl_t* buffer = new hvl_t[dims[1]];
    hid_t memtype = H5Tvlen_create(H5T_NATIVE_FLOAT);

    hsize_t dataspace_start[3] = {0, 0};
    hsize_t dataspace_count[3] = {1, dims[1]};

    std::vector<boost::multi_array<float, 2> > mu(spliced_tgroup_indexes.size());
    for (size_t i = 0; i < dims[1]; ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        mu[i].resize(boost::extents[dims[0]][tgroup_tids[tgroup].size()]);
    }

    for (size_t i = 0; i < num_samples; ++i) {
        dataspace_start[0] = i;

        H5Sselect_hyperslab_checked(dataspace, H5S_SELECT_SET,
                                    dataspace_start, NULL,
                                    dataspace_count, NULL);

        H5Dread_checked(dataset, memtype, mem_dataspace, dataspace,
                        H5P_DEFAULT, buffer);
        float* xs = reinterpret_cast<float*>(buffer[stg].p);
        output[i] = xs[idx];
    }

    return output;
}


std::vector<std::vector<float> > Summarize::transcript_condition_splicing(const char* transcript_id)
{
    std::vector<std::vector<float> > output(C);
    for (size_t i = 0; i < C; ++i) {
        output[i].resize(num_samples);
    }

    size_t tid = 0;
    for (; tid < transcript_ids.size(); ++tid) {
        if (transcript_ids[tid] == transcript_id) break;
    }

    if (tid == transcript_ids.size()) {
        Logger::abort("Transcript ID %s not present in dataset.", transcript_id);
    }

    size_t tg = tgroup[tid];

    size_t stg = 0;
    for (; stg < spliced_tgroup_indexes.size(); ++stg) {
        if (spliced_tgroup_indexes[stg] == tg) break;
    }

    size_t idx = 0;
    for (; idx < tgroup_tids[tg].size(); ++idx) {
        if (tgroup_tids[tg][idx] == tid) break;
    }

    typedef boost::multi_array<float, 3> marray_t;
    std::vector<marray_t> splicing(spliced_tgroup_indexes.size());
    condition_splicing(splicing);

    for (unsigned int j = 0; j < C; ++j) {
        for (size_t k = 0; k < num_samples; ++k) {
            output[j][k] = splicing[stg][k][j][idx];
        }
    }

    return output;
}


void Summarize::differential_feature_splicing(FILE* output,
                                              double credible_interval,
                                              double effect_size)
{
    if (isnan(effect_size)) {
        effect_size = 0.2;
    }
    bool print_credible_interval = !isnan(credible_interval);
    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    std::vector<Interval> feature_intervals;
    std::vector<std::vector<unsigned int> > including_tids;
    std::vector<std::vector<unsigned int> > excluding_tids;
    std::vector<GeneFeatureType> feature_types;
    read_gene_features(feature_intervals, including_tids, excluding_tids,
                       feature_types);

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

    // indexed by spliced tgroup, sample, condition, within tgroup index
    typedef boost::multi_array<float, 3> marray_t;
    std::vector<marray_t> splicing(spliced_tgroup_indexes.size());
    condition_splicing(splicing);

    marray_t tgroup_means;
    read_tgroup_mean(tgroup_means);

    std::vector<float> work(num_samples);

    fprintf(output,
            "gene_names"
            "\tgene_ids"
            "\tincluding_transcript_ids"
            "\texcluding_transcript_ids"
            "\tlocus"
            "\ttype"
            "\tcondition_a"
            "\tcondition_b"
            "\tdown_pr"
            "\tup_pr"
            "\tmedian_change");

    if (print_credible_interval) {
        fprintf(output, "\tlower_log2_fold_change\tupper_log2_fold_change");
    }
    fputc('\n', output);

    boost::multi_array<float, 2> spliced_in_proportion(
            boost::extents[C][num_samples]);

    std::set<GeneName> feature_gene_names;
    std::set<GeneID> feature_gene_ids;
    std::set<TranscriptID> feature_including_tids, feature_excluding_tids;
    std::set<unsigned int> used_tgroups, including_tgroups, excluding_tgroups;
    bool first_item;

    for (size_t i = 0; i < feature_intervals.size(); ++i) {
        // whitelist tgroups that actually both spliceforms
        including_tgroups.clear();
        excluding_tgroups.clear();
        used_tgroups.clear();
        BOOST_FOREACH (unsigned int tid, including_tids[i]) {
            including_tgroups.insert(tgroup[tid]);
        }

        BOOST_FOREACH (unsigned int tid, excluding_tids[i]) {
            excluding_tgroups.insert(tgroup[tid]);
        }

        std::set_intersection(including_tgroups.begin(), including_tgroups.end(),
                              excluding_tgroups.begin(), excluding_tgroups.end(),
                              std::inserter(used_tgroups, used_tgroups.begin()));

        if (used_tgroups.empty()) continue;

        // pre-compute spliced in proportions
        for (unsigned int j = 0; j < C; ++j) {
            for (size_t k = 0; k < num_samples; ++k) {
                double spliced_in_tgroup_expr_total = 0.0;
                double spliced_in = -INFINITY;
                BOOST_FOREACH (unsigned int tid, including_tids[i]) {
                    unsigned int tg = tgroup[tid];
                    if (used_tgroups.find(tg) == used_tgroups.end()) {
                        continue;
                    }

                    unsigned int stg = tgroup_spliced_tgroup[tg];
                    if (stg == (unsigned int) -1) {
                        continue;
                    }
                    unsigned int stid = tid_tgroup_index[tid];

                    double tgm = tgroup_means[k][j][tg];
                    if (spliced_in_tgroup_expr_total == 0.0) {
                        spliced_in_tgroup_expr_total = tgm;
                    }
                    else {
                        spliced_in_tgroup_expr_total = logaddexp(
                                spliced_in_tgroup_expr_total, tgm);
                    }

                    double x = log(splicing[stg][k][j][stid]) + tgm;

                    if (isinf(spliced_in)) spliced_in = x;
                    else spliced_in = logaddexp(spliced_in, x);
                }
                spliced_in -= spliced_in_tgroup_expr_total;

                double spliced_out = -INFINITY;
                double spliced_out_tgroup_expr_total = 0.0;
                BOOST_FOREACH (unsigned int tid, excluding_tids[i]) {
                    unsigned int tg = tgroup[tid];
                    if (used_tgroups.find(tg) == used_tgroups.end()) continue;

                    unsigned int stg = tgroup_spliced_tgroup[tg];
                    if (stg == (unsigned int) -1) continue;
                    unsigned int stid = tid_tgroup_index[tid];

                    double tgm = tgroup_means[k][j][tg];
                    if (spliced_out_tgroup_expr_total == 0.0) {
                        spliced_out_tgroup_expr_total = tgm;
                    }
                    else {
                        spliced_out_tgroup_expr_total = logaddexp(
                                spliced_out_tgroup_expr_total, tgm);
                    }

                    double x = log(splicing[stg][k][j][stid]) +
                               tgroup_means[k][j][tg];
                    if (isinf(spliced_out)) spliced_out = x;
                    else spliced_out = logaddexp(spliced_out, x);
                }
                spliced_out -= spliced_out_tgroup_expr_total;

                spliced_in_proportion[j][k] =
                    exp(spliced_in - logaddexp(spliced_in, spliced_out));
            }
        }

        feature_gene_names.clear();
        feature_gene_ids.clear();
        feature_including_tids.clear();
        feature_excluding_tids.clear();

        BOOST_FOREACH (unsigned int tid, including_tids[i]) {
            if (used_tgroups.find(tgroup[tid]) == used_tgroups.end()) {
                continue;
            }
            feature_gene_names.insert(gene_names[tid]);
            feature_gene_ids.insert(gene_ids[tid]);
            feature_including_tids.insert(transcript_ids[tid]);
        }

        BOOST_FOREACH (unsigned int tid, excluding_tids[i]) {
            if (used_tgroups.find(tgroup[tid]) == used_tgroups.end()) {
                continue;
            }
            feature_gene_names.insert(gene_names[tid]);
            feature_gene_ids.insert(gene_ids[tid]);
            feature_excluding_tids.insert(transcript_ids[tid]);
        }

        for (unsigned int condition_a = 0; condition_a < C - 1; ++condition_a) {
            for (unsigned int condition_b = condition_a + 1; condition_b < C; ++condition_b) {
                // gene_names
                first_item = true;
                BOOST_FOREACH (const GeneName& gene_name, feature_gene_names) {
                    if (!first_item) fputc(',', output);
                    else first_item = false;
                    fputs(gene_name.get().c_str(), output);
                }

                // gene_ids
                fputc('\t', output);
                first_item = true;
                BOOST_FOREACH (const GeneID& gene_id, feature_gene_ids) {
                    if (!first_item) fputc(',', output);
                    else first_item = false;
                    fputs(gene_id.get().c_str(), output);
                }

                // transcript_ids
                fputc('\t', output);
                first_item = true;
                BOOST_FOREACH (const TranscriptID& transcript_id, feature_including_tids) {
                    if (!first_item) fputc(',', output);
                    else first_item = false;
                    fputs(transcript_id.get().c_str(), output);
                }

                fputc('\t', output);
                first_item = true;
                BOOST_FOREACH (const TranscriptID& transcript_id, feature_excluding_tids) {
                    if (!first_item) fputc(',', output);
                    else first_item = false;
                    fputs(transcript_id.get().c_str(), output);
                }

                // locus
                fputc('\t', output);
                fprintf(output, "%s:%ld-%ld(%c)",
                        feature_intervals[i].seqname.get().c_str(),
                        feature_intervals[i].start,
                        feature_intervals[i].end,
                        feature_intervals[i].strand == strand_pos ? '+' :
                        feature_intervals[i].strand == strand_neg ? '-' : '.');

                // type
                fputc('\t', output);
                if (feature_types[i] == GENE_FEATURE_CASSETTE_EXON) {
                    fputs("cassette_exon", output);
                }
                else if (feature_types[i] == GENE_FEATURE_RETAINED_INTRON) {
                    fputs("retained_intron", output);
                }
                else {
                    fputs("unknown_feature", output);
                }

                // condition_a, condition_b
                fprintf(output, "\t%s\t%s",
                        condition_names[condition_a].c_str(),
                        condition_names[condition_b].c_str());

                double down_pr = 0.0, up_pr = 0.0;
                std::fill(work.begin(), work.end(), 0.0);
                for (size_t k = 0; k < num_samples; ++k) {
                    work[k] = spliced_in_proportion[condition_a][k] -
                              spliced_in_proportion[condition_b][k];
                    if (fabs(work[k]) >= effect_size) {
                        if (work[k] <= 0) down_pr += 1;
                        else up_pr += 1;
                    }
                }
                down_pr /= num_samples;
                up_pr /= num_samples;

                std::sort(work.begin(), work.end());
                fprintf(output, "\t%0.3f\t%0.3f\t%e",
                        down_pr, up_pr, work[num_samples/2]);
                if (print_credible_interval) {
                    fprintf(output, "\t%e\t%e",
                            work[lround((num_samples - 1) * lower_quantile)],
                            work[lround((num_samples - 1) * upper_quantile)]);
                }
                fputc('\n', output);
            }
        }
    }
}


void Summarize::differential_gene_expression(FILE* output, double credible_interval,
                                             double effect_size)
{
    if (isnan(effect_size)) {
        effect_size = 1.0;
    }
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    fprintf(output,
            "gene_name\tgene_id\ttranscript_ids\tcondition_a\tcondition_b\t"
            "down_pr\tup_pr\tmedian_log2_fold_change");
    if (print_credible_interval) {
        fprintf(output, "\tlower_log2_fold_change\tupper_log2_fold_change");
    }
    fputc('\n', output);


    typedef std::pair<GeneID, GeneName> item_t;
    size_t num_genes = gid_to_tids.size();
    boost::multi_array<float, 3> expr_data(boost::extents[num_samples][C][num_genes]);
    condition_gene_expression(expr_data);
    std::vector<double> work(num_samples);

    for (unsigned int condition_a = 0; condition_a < C - 1; ++condition_a) {
        for (unsigned int condition_b = condition_a + 1; condition_b < C; ++condition_b) {
            size_t i = 0;
            BOOST_FOREACH (item_t item, gid_to_gene_name) {
                GeneID gene_id = item.first;
                double down_pr = 0, up_pr = 0;

                for (unsigned int j = 0; j < num_samples; ++j) {
                    work[j] = log2(expr_data[j][condition_a][i] /
                                   expr_data[j][condition_b][i]);
                    if (work[j] >= effect_size) {
                        up_pr += 1;
                    }
                    else if (work[j] <= -effect_size) {
                        down_pr += 1;
                    }
                }
                up_pr /= num_samples;
                down_pr /= num_samples;
                std::sort(work.begin(), work.end());

                fprintf(output, "%s\t%s\t",
                        gid_to_gene_name[gene_id].get().c_str(),
                        gene_id.get().c_str());

                unsigned int l = 0;
                BOOST_FOREACH (unsigned int tid, gid_to_tids[gene_id]) {
                    if (l > 0) {
                        fputc(',', output);
                    }
                    fputs(transcript_ids[tid].get().c_str(), output);
                    ++l;
                }

                fprintf(output, "\t%s\t%s\t%0.3f\t%0.3f\t%e",
                        condition_names[condition_a].c_str(),
                        condition_names[condition_b].c_str(),
                        down_pr, up_pr, work[num_samples/2]);

                if (print_credible_interval) {
                    fprintf(output, "\t%e\t%e",
                            work[lround((num_samples - 1) * lower_quantile)],
                            work[lround((num_samples - 1) * upper_quantile)]);
                }
                fputc('\n', output);

                ++i;
            }
        }
    }
}


void Summarize::differential_transcript_expression(FILE* output, double credible_interval,
                                                   double effect_size)
{
    if (isnan(effect_size)) {
        effect_size = 1.0;
    }
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    fprintf(output,
            "gene_name\tgene_id\ttranscript_id\tcondition_a\tcondition_b\t"
            "down_pr\tup_pr\tmedian_log2_fold_change");
    if (print_credible_interval) {
        fprintf(output, "\tlower_log2_fold_change\tupper_log2_fold_change");
    }
    fputc('\n', output);

    boost::multi_array<float, 3> expr_data(boost::extents[num_samples][C][N]);
    condition_transcript_expression(expr_data);
    std::vector<double> work(num_samples);

    for (unsigned int condition_a = 0; condition_a < C - 1; ++condition_a) {
        for (unsigned int condition_b = condition_a + 1; condition_b < C; ++condition_b) {
            for (unsigned int i = 0; i < N; ++i) {
                double down_pr = 0, up_pr = 0;
                for (unsigned int j = 0; j < num_samples; ++j) {
                    work[j] = log2(expr_data[j][condition_a][i] /
                                   expr_data[j][condition_b][i]);
                    if (work[j] >= effect_size) {
                        up_pr += 1;
                    }
                    else if (work[j] <= -effect_size) {
                        down_pr += 1;
                    }
                }
                up_pr /= num_samples;
                down_pr /= num_samples;
                std::sort(work.begin(), work.end());

                fprintf(output, "%s\t%s\t%s\t%s\t%s\t%0.3f\t%0.3f\t%e",
                        gene_names[i].get().c_str(),
                        gene_ids[i].get().c_str(),
                        transcript_ids[i].get().c_str(),
                        condition_names[condition_a].c_str(),
                        condition_names[condition_b].c_str(),
                        down_pr, up_pr, work[num_samples/2]);
                if (print_credible_interval) {
                    fprintf(output, "\t%e\t%e",
                            work[lround((num_samples - 1) * lower_quantile)],
                            work[lround((num_samples - 1) * upper_quantile)]);
                }
                fputc('\n', output);
            }
        }
    }
}


void Summarize::condition_transcript_expression(FILE* output,
                                                double credible_interval)
{
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    fprintf(output, "gene_name\tgene_id\ttranscript_id");
    for (unsigned int i = 0; i < C; ++i) {
        fprintf(output, "\t%s_adjusted_tpm",
                metadata.sample_conditions[i].c_str());
        if (print_credible_interval) {
            fprintf(output,
                    "\t%s_adjusted_tpm_lower\t%s_adjusted_tpm_upper",
                    metadata.sample_conditions[i].c_str(),
                    metadata.sample_conditions[i].c_str());
        }
    }
    fputc('\n', output);

    boost::multi_array<float, 3> expr_data(boost::extents[num_samples][C][N]);
    condition_transcript_expression(expr_data);
    std::vector<double> work(num_samples);

    for (unsigned int i = 0; i < N; ++i) {
        fprintf(output, "%s\t%s\t%s",
                gene_names[i].get().c_str(),
                gene_ids[i].get().c_str(),
                transcript_ids[i].get().c_str());

        for (unsigned int j = 0; j < C; ++j) {
            for (unsigned int k = 0; k < num_samples; ++k) {
                work[k] = expr_data[k][j][i];
            }
            std::sort(work.begin(), work.end());
            fprintf(output, "\t%e", 1e6 * work[num_samples / 2]);
            if (print_credible_interval) {
                fprintf(output, "\t%e\t%e",
                        1e6 * work[lround((num_samples - 1) * lower_quantile)],
                        1e6 * work[lround((num_samples - 1) * upper_quantile)]);
            }
        }
        fputc('\n', output);
    }
}


// output indexed by sample number, condition, tid
void Summarize::condition_transcript_expression(boost::multi_array<float, 3>& output)
{
    // indexed by: spliced tgroup, sample number, condition, within tgroup tid
    std::vector<boost::multi_array<float, 3> > splicing_data(spliced_tgroup_indexes.size());
    condition_splicing(splicing_data);

    // indexed by: sample number, condition, tgroup
    boost::multi_array<float, 3> tgroup_mean_data;
    read_tgroup_mean(tgroup_mean_data);

    tgroup_mean_data.begin();

    for (size_t i = 0; i < num_samples; ++i) {
        for (size_t j = 0; j < C; ++j) {
            for (size_t k = 0; k < T; ++k) {
                tgroup_mean_data[i][j][k] = exp(tgroup_mean_data[i][j][k]);
            }
        }
    }

    // fill in alternatively spliced transcripts
    for (unsigned int stg = 0; stg < spliced_tgroup_indexes.size(); ++stg) {
        unsigned int tg = spliced_tgroup_indexes[stg];
        for (unsigned int k = 0; k < tgroup_tids[tg].size(); ++k) {
            for (size_t i = 0; i < num_samples; ++i) {
                for (size_t j = 0; j < C; ++j) {
                    output[i][j][k] =
                        tgroup_mean_data[i][j][tg] *
                        splicing_data[stg][i][j][k];
                }
            }
        }
    }

    // fill in transcripts that are not alternatively spliced
    for (size_t k = 0; k < N; ++k) {
        unsigned int tg = tgroup[k];
        if (tgroup_tids[tg].size() != 1) continue;

        for (size_t i = 0; i < num_samples; ++i) {
            for (size_t j = 0; j < C; ++j) {
                output[i][j][tgroup_tids[tg][k]] = tgroup_mean_data[i][j][tg];
            }
        }
    }
}


void Summarize::condition_gene_expression(FILE* output,
                                          double credible_interval)
{
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    fprintf(output, "gene_name\tgene_id\ttranscript_ids");
    for (unsigned int i = 0; i < C; ++i) {
        fprintf(output, "\t%s_adjusted_tpm",
                metadata.sample_conditions[i].c_str());
        if (print_credible_interval) {
            fprintf(output,
                    "\t%s_adjusted_tpm_lower\t%s_adjusted_tpm_upper",
                    metadata.sample_conditions[i].c_str(),
                    metadata.sample_conditions[i].c_str());
        }
    }
    fputc('\n', output);

    size_t num_genes = gid_to_tids.size();
    boost::multi_array<float, 3> expr_data(boost::extents[num_samples][C][num_genes]);
    condition_gene_expression(expr_data);
    std::vector<double> work(num_samples);

    unsigned int i = 0;
    typedef std::pair<GeneID, GeneName> item_t;
    BOOST_FOREACH (item_t item, gid_to_gene_name) {
        GeneID gene_id = item.first;
        fprintf(output, "%s\t%s\t",
                gid_to_gene_name[gene_id].get().c_str(),
                gene_id.get().c_str());

        unsigned int l = 0;
        BOOST_FOREACH (unsigned int tid, gid_to_tids[gene_id]) {
            if (l > 0) {
                fputc(',', output);
            }
            fputs(transcript_ids[tid].get().c_str(), output);
            ++l;
        }

        for (unsigned int j = 0; j < C; ++j) {
            for (unsigned int k = 0; k < num_samples; ++k) {
                work[k] = expr_data[k][j][i];
            }
            std::sort(work.begin(), work.end());
            fprintf(output, "\t%e", 1e6 * work[num_samples / 2]);
            if (print_credible_interval) {
                fprintf(output, "\t%e\t%e",
                        1e6 * work[lround((num_samples - 1) * lower_quantile)],
                        1e6 * work[lround((num_samples - 1) * upper_quantile)]);
            }
        }
        fputc('\n', output);
        ++i;
    }
}


// output index by sample number, condition, gid
void Summarize::condition_gene_expression(boost::multi_array<float, 3>& output)
{
    boost::multi_array<float, 3> transcript_expr_data(boost::extents[num_samples][C][N]);
    condition_transcript_expression(transcript_expr_data);

    typedef std::pair<GeneID, GeneName> item_t;
    for (size_t i = 0; i < num_samples; ++i) {
        for (size_t j = 0; j < C; ++j) {
            size_t k = 0;
            BOOST_FOREACH (item_t item, gid_to_gene_name) {
                GeneID gene_id = item.first;
                output[i][j][k] = 0.0;
                BOOST_FOREACH (size_t l, gid_to_tids[gene_id]) {
                    output[i][j][k] += transcript_expr_data[i][j][l];
                }
                ++k;
            }
        }
    }
}


void Summarize::condition_splicing(FILE* output, double credible_interval)
{
    bool print_credible_interval = !isnan(credible_interval);

    double lower_quantile = 0.5 - credible_interval/2;
    double upper_quantile = 0.5 + credible_interval/2;

    // indexed by sample_num, condition, and within-tgroup transcript
    typedef boost::multi_array<float, 3> marray_t;

    // indexed by spliced tgroup
    std::vector<marray_t> splicing(spliced_tgroup_indexes.size());

    condition_splicing(splicing);

    fprintf(output, "gene_names\tgene_ids\ttranscription_group\ttranscript_id");
    for (size_t i = 0; i < C; ++i) {
        if (print_credible_interval) {
            fprintf(output, "\t%s_splice_rate\t%s_splice_rate_lower\t%s_splice_rate_upper",
                    condition_names[i].c_str(), condition_names[i].c_str(),
                    condition_names[i].c_str());
        }
        else {
            fprintf(output, "\t%s_splice_rate", condition_names[i].c_str());
        }
    }
    fputc('\n', output);

    std::set<GeneID> tgroup_gene_ids;
    std::set<GeneName> tgroup_gene_names;

    std::vector<float> work(num_samples);
    for (size_t i = 0; i < spliced_tgroup_indexes.size(); ++i) {
        size_t tgroup = spliced_tgroup_indexes[i];
        tgroup_gene_names.clear();
        tgroup_gene_ids.clear();
        BOOST_FOREACH (unsigned int tid, tgroup_tids[tgroup]) {
            tgroup_gene_names.insert(gene_names[tid]);
            tgroup_gene_ids.insert(gene_ids[tid]);
        }

        for (size_t j = 0; j < tgroup_tids[tgroup].size(); ++j) {
            unsigned int tid = tgroup_tids[tgroup][j];

            // gene_names
            bool first_item = true;
            BOOST_FOREACH (const GeneName& gene_name, tgroup_gene_names) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", gene_name.get().c_str());
            }
            fputc('\t', output);

            // gene_ids
            first_item = true;
            BOOST_FOREACH (const GeneID& gene_id, tgroup_gene_ids) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", gene_id.get().c_str());
            }
            fputc('\t', output);

            // transcription_group
            first_item = true;
            BOOST_FOREACH (unsigned int tid, tgroup_tids[tgroup]) {
                if (!first_item) fputc(',', output);
                else first_item = false;
                fprintf(output, "%s", transcript_ids[tid].get().c_str());
            }
            fputc('\t', output);

            // transcript_id
            fprintf(output, "%s", transcript_ids[tid].get().c_str());

            for (size_t k = 0; k < C; ++k) {
                for (size_t l = 0; l < num_samples; ++l) {
                    work[l] = splicing[i][l][k][j];
                }
                std::sort(work.begin(), work.end());

                if (print_credible_interval) {
                    fprintf(output, "\t%f\t%f\t%f",
                            (double) work[lround((num_samples -1) * lower_quantile)],
                            (double) work[num_samples/2],
                            (double) work[lround((num_samples -1) * upper_quantile)]);
                }
                else {
                    fprintf(output, "\t%f", (double) work[num_samples/2]);
                }
            }
            fputc('\n', output);
        }
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

    hid_t dataset = H5Dopen2_checked(h5_file, "/condition/tgroup_mean", H5P_DEFAULT);

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

        H5Dread_checked(dataset, H5T_NATIVE_FLOAT,
                        mem_dataspace, dataspace, H5P_DEFAULT,
                        &condition_data_a[j].at(0));
    }

    // read ccondition_b data
    for (size_t j = 0; j < num_samples; ++j) {
        file_dataspace_start[0] = j;
        file_dataspace_start[1] = condition_b;

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            file_dataspace_start, NULL,
                            file_dataspace_dims, NULL);

        H5Dread_checked(dataset, H5T_NATIVE_FLOAT,
                        mem_dataspace, dataspace, H5P_DEFAULT,
                        &condition_data_b[j].at(0));
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
                            transcript_ids[tid].get().c_str(), gene_ids[tid].get().c_str(),
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
    hid_t dataset = H5Dopen2_checked(h5_file, "/transcript_quantification", H5P_DEFAULT);

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
    H5Dread_checked(dataset, H5T_NATIVE_FLOAT, dataspace, dataspace,
                    H5P_DEFAULT, &Q.data()[0]);

    fprintf(output, "transcript_id");
    for (size_t i = 0; i < K; ++i) {
        fprintf(output, "\tsample%lu", (unsigned long) i);
    }
    fputc('\n', output);

    for (size_t i = 0; i < N; ++i) {
        fputs(transcript_ids[i].get().c_str(), output);
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


void Summarize::read_gene_features(std::vector<Interval>& feature_intervals,
                                    std::vector<std::vector<unsigned int> >& including_tids,
                                    std::vector<std::vector<unsigned int> >& excluding_tids,
                                    std::vector<GeneFeatureType>& feature_types)
{
    hid_t dataset = H5Dopen2_checked(h5_file, "/features/seqname", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    size_t n = dims[0];

    // read sequence names
    hid_t varstring_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(varstring_type, H5T_VARIABLE);
    std::vector<char*> seqnames(n);
    H5Dread_checked(dataset, varstring_type, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, &seqnames.at(0));
    H5Dclose(dataset);

    // read starts
    dataset = H5Dopen2_checked(h5_file, "/features/start", H5P_DEFAULT);
    std::vector<long> starts(n);
    H5Dread_checked(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &starts.at(0));
    H5Dclose(dataset);

    // read ends
    dataset = H5Dopen2_checked(h5_file, "/features/end", H5P_DEFAULT);
    std::vector<long> ends(n);
    H5Dread_checked(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &ends.at(0));
    H5Dclose(dataset);

    // read strands
    dataset = H5Dopen2_checked(h5_file, "/features/strand", H5P_DEFAULT);
    std::vector<char> strands(n);
    H5Dread_checked(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &strands.at(0));
    H5Dclose(dataset);

    // read including_tids
    dataset = H5Dopen2_checked(h5_file, "/features/including_tid", H5P_DEFAULT);
    hid_t tids_type = H5Tvlen_create(H5T_NATIVE_UINT);
    std::vector<hvl_t> including_tids_data(n);
    H5Dread_checked(dataset, tids_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &including_tids_data.at(0));
    H5Dclose(dataset);

    // read excluding_tids
    dataset = H5Dopen2_checked(h5_file, "/features/excluding_tid", H5P_DEFAULT);
    std::vector<hvl_t> excluding_tids_data(n);
    H5Dread_checked(dataset, tids_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &excluding_tids_data.at(0));
    H5Dclose(dataset);

    // read types
    dataset = H5Dopen2_checked(h5_file, "/features/type", H5P_DEFAULT);
    std::vector<unsigned int> types_data(n);
    H5Dread_checked(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &types_data.at(0));
    H5Dclose(dataset);

    // construct the data
    feature_intervals.resize(n);
    including_tids.resize(n);
    excluding_tids.resize(n);
    feature_types.resize(n);

    for (size_t i = 0; i < n; ++i) {
        feature_intervals[i].seqname = std::string(seqnames[i]);
        feature_intervals[i].start = starts[i];
        feature_intervals[i].end = ends[i];
        if (strands[i] == '+') {
            feature_intervals[i].strand = strand_pos;
        }
        else if (strands[i] == '-') {
            feature_intervals[i].strand = strand_neg;
        }
        else {
            feature_intervals[i].strand = strand_na;
        }

        feature_types[i] = (GeneFeatureType) types_data[i];

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
    hid_t dataset = H5Dopen2_checked(h5_file, "/condition/tgroup_mean", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[3]; // dims are: num_samples, C, T
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    output.resize(boost::extents[dims[0]][dims[1]][dims[2]]);

    H5Dread_checked(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, output.data());
    H5Dclose(dataset);
}


Summarize::~Summarize()
{
    H5Fclose(h5_file);
}


void Summarize::read_metadata(IsolatorMetadata& metadata)
{
    hid_t group = H5Gopen2(h5_file, "/metadata", H5P_DEFAULT);
    if (group < 0) {
        Logger::abort("Failed to open the '/metadata' group.");
    }

    hid_t varstring_type = H5Tcopy(H5T_C_S1);
    if (varstring_type < 0 || H5Tset_size(varstring_type, H5T_VARIABLE) < 0) {
        Logger::abort("HDF5 type creation failed.");
    }

    hid_t attr;
    std::vector<char*> data(1);

    attr = H5Aopen_checked(group, "command_line", H5P_DEFAULT);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.command_line = data[0];
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "version", H5P_DEFAULT);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.version = data[0];
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "commit", H5P_DEFAULT);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.commit = data[0];
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "date", H5P_DEFAULT);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.date = data[0];
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "elapsed_seconds", H5P_DEFAULT);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.elapsed_seconds = data[0];
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "sample_filenames", H5P_DEFAULT);
    hid_t samples_dataspace = H5Aget_space(attr);
    hsize_t num_samples = 0;
    H5Sget_simple_extent_dims(samples_dataspace, &num_samples, NULL);
    H5Sclose(samples_dataspace);
    data.resize(num_samples);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.sample_filenames.resize(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        metadata.sample_filenames[i] = data[i];
    }
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "sample_names", H5P_DEFAULT);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.sample_names.resize(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        metadata.sample_names[i] = data[i];
    }
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "sample_conditions", H5P_DEFAULT);
    H5Aread(attr, varstring_type, &data.at(0));
    metadata.sample_conditions.resize(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        metadata.sample_conditions[i] = data[i];
    }
    H5Aclose(attr);

    attr = H5Aopen_checked(group, "rng_seed", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &metadata.rng_seed);
    H5Aclose(attr);

    H5Tclose(varstring_type);
    H5Gclose(group);
}


const std::vector<TranscriptID>& Summarize::get_transcript_ids()
{
    return transcript_ids;
}


const std::vector<GeneID>& Summarize::get_gene_ids()
{
    return gene_ids;
}


const std::vector<unsigned int>& Summarize::get_tgroups()
{
    return tgroup;
}



