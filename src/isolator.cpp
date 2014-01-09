
//#include <google/profiler.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <getopt.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <map>
#include <sys/time.h>
#include <unistd.h>

#include "analyze.hpp"
#include "config.h"
#include "constants.hpp"
#include "fragment_model.hpp"
#include "linalg.hpp"
#include "logger.hpp"
#include "sampler.hpp"
#include "summarize.hpp"
#include "transcripts.hpp"

using namespace boost::accumulators;

const char* isolator_logo =
"     _           _       _\n"
"    (_)___  ___ | | __ _| |_ ___  _ __\n"
"    | / __|/ _ \\| |/ _` | __/ _ \\| '__|\n"
"    | \\__ \\ (_) | | (_| | || (_) | |\n"
"    |_|___/\\___/|_|\\__,_|\\__\\___/|_|\n";


static void print_logo()
{
    printf("%s\n     Version: %s\n     Instruction set: %s\n\n",
           isolator_logo, VERSION, LINALG_INSTR_SET);
}


void print_describe_usage(FILE* fout)
{
    fprintf(fout, "Usage: isolator describe isolator-output.h5\n");
}


void print_describe_help(FILE* fout)
{
    print_describe_usage(fout);
    fprintf(fout,
            "\nOptions:\n"
            "-h, --help                Print this help message\n"
            "-o, --output=FILE         Write output to the given file (default: standard out)\n"
            "\n"
            "See 'isolator help summarize' for more.\n");
}


int isolator_describe(int argc, char* argv[])
{
    Logger::start();

    static struct option long_options[] =
    {
        {"help",        no_argument,       NULL, 'h'},
        {"output",      required_argument, NULL, 'o'},
        {0, 0, 0, 0}
    };

    int opt;
    int opt_idx;
    const char* output_filename = NULL;

    while (true) {
        opt = getopt_long(argc, argv, "ho:", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_describe_help(stdout);
                return 0;

            case 'o':
                output_filename = optarg;
                break;

            case '?':
                fprintf(stderr, "\n");
                print_describe_help(stderr);
                return 1;

            default:
                abort();
        }
    }

    /* no positional argumens */
    if (optind == argc) {
        print_describe_usage(stdout);
        return 0;
    }

    /* too many */
    else if (optind + 1 > argc) {
        fprintf(stderr, "Too many arguments.\n\n");
        print_describe_usage(stderr);
        return 1;
    }

    const char* input_filename = argv[optind];
    Summarize summarize(input_filename);

    FILE* output_file = stdout;
    if (output_filename) {
        output_file = fopen(output_filename, "w");
        if (!output_file) {
            Logger::abort("Can't open file %s for writing.\n", output_filename);
        }
    }

    IsolatorMetadata metadata;
    summarize.read_metadata(metadata);

    const std::vector<TranscriptID>& transcript_ids = summarize.get_transcript_ids();
    const std::vector<GeneID>& gene_ids = summarize.get_gene_ids();
    const std::vector<unsigned int>& tgroups = summarize.get_tgroups();

    std::set<GeneID> unique_gene_ids;
    unique_gene_ids.insert(gene_ids.begin(), gene_ids.end());

    std::set<unsigned int> unique_tgroups;
    unique_tgroups.insert(tgroups.begin(), tgroups.end());

    std::map<std::string, std::vector<std::string> > condition_samples;
    for (size_t i = 0; i < metadata.sample_filenames.size(); ++i) {
        condition_samples[metadata.sample_conditions[i]].push_back(
                metadata.sample_filenames[i]);
    }

    fprintf(output_file,
        "Experiment Design\n"
        "-----------------\n"
        "\n"
        "%lu samples in %lu conditions.\n",
        (unsigned long) metadata.sample_conditions.size(),
        (unsigned long) condition_samples.size());

    typedef std::map<std::string, std::vector<std::string> >::value_type
        condition_samples_value_type;
    BOOST_FOREACH (condition_samples_value_type& condition, condition_samples) {
        fprintf(output_file, "\nCondition %s:\n", condition.first.c_str());
        BOOST_FOREACH (std::string& filename, condition.second) {
            fprintf(output_file, "  %s\n", filename.c_str());
        }
    }

    fprintf(output_file,
        "\n\nAnalysis\n"
        "--------\n"
        "\n"
        "Gene annotations contained %lu transcripts in %lu genes and\n"
        "%lu transcription groups.\n"
        "\n"
        "Isolator version %s was run on %s taking %s seconds.\n"
        "\n"
        "Command line:\n  %s\n\n",
        (unsigned long) transcript_ids.size(),
        (unsigned long) unique_gene_ids.size(),
        (unsigned long) unique_tgroups.size(),
        metadata.version.c_str(),
        metadata.date.c_str(),
        metadata.elapsed_seconds.c_str(),
        metadata.command_line.c_str());

    return EXIT_SUCCESS;
}


void print_summarize_usage(FILE* fout)
{
    fprintf(fout, "Usage: isolator summarize strategy [options] analyze_output.h5\n");
}


void print_summarize_help(FILE* fout)
{
    print_summarize_usage(fout);
    fprintf(fout,
            "\nOptions:\n"
            "-h, --help                Print this help message\n"
            "-v, --verbose             Print a bunch of information useful mainly for debugging.\n"
            "-o, --out=FILE            Output summary to the given file. (default: standard out).\n"
            "-l, --list                List summarization strategies.\n"
            "-a, --condition_a         First condition for pairwise comparisons.\n"
            "-b, --condition_b         First condition for pairwise comparisons.\n"
            "\n"
            "See 'isolator help summarize' for more.\n");
}


void print_summarize_strategies(FILE* fout)
{
    fprintf(fout,
           "Available strategies:\n"
           "  median_transcript_expression\n"
           "  median_gene_expression\n"
           "  condition_splicing\n"
           "  condition_pairwise_splicing\n"
           "  condition_tgroup_mean\n"
           "  experiment_tgroup_sd\n"
           "  tgroup_fold_change\n"
           "  expression_samples\n"
           "  cassette_exons\n");
}


int isolator_summarize(int argc, char* argv[])
{
    Logger::start();

    static struct option long_options[] =
    {
        {"help",         no_argument,       NULL, 'h'},
        {"verbose",      no_argument,       NULL, 'v'},
        {"out",          required_argument, NULL, 'o'},
        {"list",         no_argument,       NULL, 'l'},
        {"strategy",     required_argument, NULL, 's'},
        {"condition_a",  required_argument, NULL, 'a'},
        {"condition_b",  required_argument, NULL, 'b'},
        {"unnormalized", no_argument,       NULL, 'u'},
        {"credible",     required_argument, NULL, 'c'},
        {"effect-size",  required_argument, NULL, 'e'},
        {0, 0, 0, 0}
    };

    const char* out_fn = "isolator_summarize.tsv";
    Logger::level logger_level = Logger::INFO;
    unsigned int condition_a = 0, condition_b = 1;

    int opt;
    int opt_idx;
    double credible_interval = NAN;
    double minimum_effect_size = 1.5;
    bool unnormalized = false;

    while (true) {
        opt = getopt_long(argc, argv, "hvo:s:l", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_summarize_help(stdout);
                return 0;

            case 'v':
                logger_level = Logger::DEBUG;
                break;

            case 'o':
                out_fn = optarg;
                break;

            case 'l':
                print_summarize_strategies(stdout);
                return 0;

            case 'a':
                condition_a = strtoul(optarg, NULL, 10);
                break;

            case 'b':
                condition_b = strtoul(optarg, NULL, 10);
                break;

            case 'u':
                unnormalized = true;
                break;

            case 'c':
                credible_interval = atof(optarg);
                break;

            case 'e':
                minimum_effect_size = atof(optarg);
                break;

            case '?':
                fprintf(stderr, "\n");
                print_summarize_help(stderr);
                return 1;

            default:
                abort();
        }
    }

    /* too few positional argumens */
    if (optind + 1 >= argc) {
        fprintf(stderr, "Too few arguments.\n\n");
        print_summarize_usage(stdout);
        return 0;
    }

    /* too many */
    else if (optind + 2 < argc) {
        fprintf(stderr, "Too many arguments.\n\n");
        print_summarize_usage(stderr);
        return 1;
    }

    FILE* out_file = stdout;
    if (out_fn) {
        out_file = fopen(out_fn, "w");
        if (out_file == NULL) {
            Logger::abort("Can't open file %s for writing.\n", out_fn);
        }
    }

    char* strategy = argv[optind];
    const char* in_fn = argv[optind+1];

    // strategies are not case sensitive
    for (char* c = strategy; *c; ++c) {
        *c = tolower(*c);
    }

    /* TODO: Plan for summarization strategies.
     *
     *     transcript_expression
     *     gene_expression
     *     tgroup_expression
     *
     *     differential_transcription
     *     differential_splicing
     */

    minimum_effect_size = log2(minimum_effect_size);

    Summarize summarize(in_fn);

    if (strcmp(strategy, "transcript_expression") == 0) {
        summarize.median_transcript_expression(out_file, credible_interval,
                                               unnormalized);
    }
    else if (strcmp(strategy, "gene_expression") == 0) {
        summarize.median_gene_expression(out_file, credible_interval,
                                         unnormalized);
    }
    else if (strcmp(strategy, "differential_transcription") == 0) {
        summarize.differential_transcription(out_file, credible_interval,
                                             minimum_effect_size);
    }
    else if (strcmp(strategy, "differential_splicing") == 0) {
        summarize.differential_splicing(out_file, credible_interval,
                                        minimum_effect_size);
    }


#if 0
    if (strcmp(strategy, "median_transcript_expression") == 0) {
        summarize.median_transcript_expression(out_f);
    }
    else if (strcmp(strategy, "median_gene_expression") == 0) {
        summarize.median_gene_expression(out_f);
    }
    else if (strcmp(strategy, "condition_splicing") == 0) {
        summarize.condition_splicing(out_f);
    }
    else if (strcmp(strategy, "condition_pairwise_splicing") == 0) {
        summarize.condition_pairwise_splicing(out_f);
    }
    else if (strcmp(strategy, "condition_tgroup_mean") == 0) {
        summarize.median_condition_tgroup_expression(out_f);
    }
    else if (strcmp(strategy, "experiment_tgroup_sd") == 0) {
        summarize.median_experiment_tgroup_sd(out_f);
    }
    else if (strcmp(strategy, "tgroup_fold_change") == 0) {
        summarize.tgroup_fold_change(out_f, condition_a, condition_b);
    }
    else if (strcmp(strategy, "expression_samples") == 0) {
        summarize.expression_samples(out_f);
    }
    else if (strcmp(strategy, "cassette_exons") == 0) {
        summarize.cassette_exon_pairwise_splicing(out_f);
    }
#endif

    else {
        fprintf(stderr, "No such summarization strategy: %s\n\n", strategy);
        print_summarize_strategies(stderr);
        return EXIT_FAILURE;
    }

    if (out_file != stdout) fclose(out_file);

    Logger::end();
    return EXIT_SUCCESS;
}


void print_analyze_usage(FILE* fout)
{
    fprintf(fout, "Usage: isolator analyze [options] genes.gtf a1.bam[,a2.bam...] [b1.bam[,b2.bam...]]\n");
}

void print_analyze_help(FILE* fout)
{
    print_analyze_usage(fout);
    fprintf(fout,
        "\nOptions:\n"
        "-h, --help                Print this help message\n"
        "-o, --output=FILE         File to write HDF5 output to (default: isolator_output.hdf5)\n"
        "    --introns             Input consists of a BED file containing introns."
        "    --exons               Input consists of a BED file containing exons."
        "-v, --verbose             Print a bunch of information useful mainly for debugging\n"
        "-g, --genomic-seq=FILE    Correct for sequence bias, given the a the sequence\n"
        "                          against which the reads are aligned, in FAST format.\n"
        "-p, --threads=N           number of threads to use.\n"
        "    --no-gc-correction    disable post-hoc GC-content adjustments.\n"
        "-N, --num-samples         generate this many samples (default: 250)\n"
        "-B, --burnin              warmup for this many samples before collecting data (default: 100)\n\n"
        "See 'isolator help teste' for more.\n");
}


static void write_metadata(hid_t file_id, const IsolatorMetadata& metadata)
{
    hid_t dataspace = H5Screate(H5S_SCALAR);
    if (dataspace < 0) {
        Logger::abort("HDF5 dataspace creation failed.");
    }

    hid_t varstring_type = H5Tcopy(H5T_C_S1);
    if (varstring_type < 0 || H5Tset_size(varstring_type, H5T_VARIABLE) < 0) {
        Logger::abort("HDF5 type creation failed.");
    }

    hid_t group = H5Gcreate1(file_id, "/metadata", 0);
    if (group < 0) {
        Logger::abort("HDF5 group creation failed.");
    }

    hid_t attr;
    const char* attr_value;

    attr = H5Acreate1(group, "command_line", varstring_type, dataspace, H5P_DEFAULT);
    attr_value = metadata.command_line.c_str();
    H5Awrite(attr, varstring_type, &attr_value);
    H5Aclose(attr);

    attr = H5Acreate1(group, "version", varstring_type, dataspace, H5P_DEFAULT);
    attr_value = metadata.version.c_str();
    H5Awrite(attr, varstring_type, &attr_value);
    H5Aclose(attr);

    attr = H5Acreate1(group, "date", varstring_type, dataspace, H5P_DEFAULT);
    attr_value = metadata.date.c_str();
    H5Awrite(attr, varstring_type, &attr_value);
    H5Aclose(attr);

    attr = H5Acreate1(group, "elapsed_seconds", varstring_type, dataspace, H5P_DEFAULT);
    attr_value = metadata.elapsed_seconds.c_str();
    H5Awrite(attr, varstring_type, &attr_value);
    H5Aclose(attr);

    hsize_t sample_dataspace_dims[1] = {metadata.sample_filenames.size()};
    hid_t sample_dataspace = H5Screate_simple(1, sample_dataspace_dims, NULL);
    if (sample_dataspace < 0) {
        Logger::abort("HDF5 dataspace creation failed.");
    }

    attr = H5Acreate1(group, "sample_filenames", varstring_type,
                      sample_dataspace, H5P_DEFAULT);
    const char** sample_attr_value = new const char*[metadata.sample_filenames.size()];
    for (size_t i = 0; i < metadata.sample_filenames.size(); ++i) {
        sample_attr_value[i] = metadata.sample_filenames[i].c_str();
    }
    H5Awrite(attr, varstring_type, sample_attr_value);
    H5Aclose(attr);

    attr = H5Acreate1(group, "sample_conditions", varstring_type,
                      sample_dataspace, H5P_DEFAULT);
    for (size_t i = 0; i < metadata.sample_conditions.size(); ++i) {
        sample_attr_value[i] = metadata.sample_conditions[i].c_str();
    }
    H5Awrite(attr, varstring_type, sample_attr_value);
    H5Aclose(attr);

    delete [] sample_attr_value;
    H5Sclose(sample_dataspace);
    H5Gclose(group);
    H5Tclose(varstring_type);
    H5Sclose(dataspace);
}


// Write a table with cassette exon information, which should look something
// like:
//   seqname, start, end, strand, inclusion_tids, exclusion_tids
//
static void write_cassette_exon_data(hid_t file_id, TranscriptSet& ts)
{
    std::vector<Interval> cassette_exons;
    std::vector<std::vector<unsigned int> > including_tids, excluding_tids;
    ts.get_cassette_exons(cassette_exons, including_tids, excluding_tids);
    size_t n = cassette_exons.size();
    Logger::debug("%lu cassette exons", (unsigned long) n);

    herr_t status;

    // create types for each field
    hid_t varstring_type = H5Tcopy(H5T_C_S1);
    if (varstring_type < 0 || H5Tset_size(varstring_type, H5T_VARIABLE) < 0) {
        Logger::abort("HDF5 type creation failed.");
    }

    hid_t tids_type = H5Tvlen_create(H5T_NATIVE_UINT);
    if (tids_type < 0) {
        Logger::abort("HDF5 type creation failed.");
    }

    if (H5Gcreate1(file_id, "/cassette_exons", 0) < 0) {
        Logger::abort("HDF5 group creation failed.");
    }

    hsize_t dims[1] = { n };
    hid_t cassette_exon_dataspace = H5Screate_simple(1, dims, NULL);
    if (cassette_exon_dataspace < 0) {
        Logger::abort("HDF5 dataspace creation failed.");
    }

    hid_t dataset_create_property = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(dataset_create_property, H5D_CHUNKED);
    H5Pset_chunk(dataset_create_property, 1, dims);
    H5Pset_deflate(dataset_create_property, 7);

    // write sequence names
    hid_t seqname_id = H5Dcreate2(file_id, "/cassette_exons/seqname",
                                  varstring_type, cassette_exon_dataspace,
                                  H5P_DEFAULT, dataset_create_property,
                                  H5P_DEFAULT);
    if (seqname_id < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    char** string_data = new char* [n];
    for (size_t i = 0; i < n; ++i) {
        const std::string& seqname = cassette_exons[i].seqname.get();
        string_data[i] = new char[seqname.size() + 1];
        memcpy(string_data[i], seqname.c_str(), seqname.size() + 1);
    }

    status = H5Dwrite(seqname_id, varstring_type, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, string_data);
    if (status < 0) Logger::abort("HDF5 write failed.");

    for (size_t i = 0; i < n; ++i) {
        delete [] string_data[i];
    }
    delete [] string_data;

    H5Dclose(seqname_id);

    // write start positions
    hid_t start_id = H5Dcreate2(file_id, "/cassette_exons/start",
                                H5T_NATIVE_LONG, cassette_exon_dataspace,
                                H5P_DEFAULT, dataset_create_property,
                                H5P_DEFAULT);

    if (start_id < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    long* long_data = new long [n];
    for (size_t i = 0; i < n; ++i) {
        long_data[i] = cassette_exons[i].start;
    }

    status = H5Dwrite(start_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, long_data);
    if (status < 0) Logger::abort("HDF5 write failed.");

    H5Dclose(start_id);

    // write end positions
    hid_t end_id = H5Dcreate2(file_id, "/cassette_exons/end",
                              H5T_NATIVE_LONG, cassette_exon_dataspace,
                              H5P_DEFAULT, dataset_create_property,
                              H5P_DEFAULT);

    if (end_id < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    for (size_t i = 0; i < n; ++i) {
        long_data[i] = cassette_exons[i].end;
    }

    status = H5Dwrite(end_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, long_data);
    if (status < 0) Logger::abort("HDF5 write failed.");

    H5Dclose(end_id);

    delete[] long_data;

    // write strands
    hid_t strand_id = H5Dcreate2(file_id, "/cassette_exons/strand",
                                 H5T_NATIVE_CHAR, cassette_exon_dataspace,
                                 H5P_DEFAULT, dataset_create_property,
                                 H5P_DEFAULT);
    if (strand_id < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    char* char_data = new char[n];
    for (size_t i = 0; i < n; ++i) {
        if (cassette_exons[i].strand == strand_pos) {
            char_data[i] = '+';
        }
        else if (cassette_exons[i].strand == strand_neg) {
            char_data[i] = '-';
        }
        else {
            char_data[i] = '.';
        }
    }

    status = H5Dwrite(strand_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, char_data);
    if (status < 0) Logger::abort("HDF5 write failed.");

    delete [] char_data;
    H5Dclose(strand_id);


    // write including tids
    hid_t including_tids_id = H5Dcreate2(file_id, "/cassette_exons/including_tids",
                                         tids_type, cassette_exon_dataspace,
                                         H5P_DEFAULT, dataset_create_property,
                                         H5P_DEFAULT);
    if (including_tids_id < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    hvl_t* vlen_data = new hvl_t [n];
    for (size_t i = 0; i < n; ++i) {
        vlen_data[i].len = including_tids[i].size();
        unsigned int* ptids = new unsigned int [including_tids[i].size()];
        for (size_t j = 0; j < including_tids[i].size(); ++j) {
            ptids[j] = including_tids[i][j];
        }
        vlen_data[i].p = reinterpret_cast<void*>(ptids);
    }

    status = H5Dwrite(including_tids_id, tids_type, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, vlen_data);
    if (status < 0) Logger::abort("HDF5 write failed.");

    H5Dclose(including_tids_id);

    // write excluding tids
    hid_t excluding_tids_id = H5Dcreate2(file_id, "/cassette_exons/excluding_tids",
                                         tids_type, cassette_exon_dataspace,
                                         H5P_DEFAULT, dataset_create_property,
                                         H5P_DEFAULT);
    if (excluding_tids_id < 0) {
        Logger::abort("HDF5 dataset creation failed.");
    }

    for (size_t i = 0; i < n; ++i) {
        vlen_data[i].len = excluding_tids[i].size();
        unsigned int* ptids = new unsigned int [excluding_tids[i].size()];
        for (size_t j = 0; j < excluding_tids[i].size(); ++j) {
            ptids[j] = excluding_tids[i][j];
        }
        delete [] reinterpret_cast<unsigned int*>(vlen_data[i].p);
        vlen_data[i].p = reinterpret_cast<void*>(ptids);
    }

    status = H5Dwrite(excluding_tids_id, tids_type, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, vlen_data);
    if (status < 0) Logger::abort("HDF5 write failed.");

    H5Dclose(excluding_tids_id);

    for (size_t i = 0; i < n; ++i) {
        delete [] reinterpret_cast<unsigned int*>(vlen_data[i].p);
    }
    delete [] vlen_data;

    H5Tclose(varstring_type);
    H5Tclose(tids_type);
}


int isolator_analyze(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"help",                 no_argument,       NULL, 'h'},
        {"output",               required_argument, NULL, 'o'},
        {"introns",              no_argument,       NULL, 0},
        {"exons",                no_argument,       NULL, 0},
        {"verbose",              no_argument,       NULL, 'v'},
        {"genomic-seq",          required_argument, NULL, 'g'},
        {"threads",              required_argument, NULL, 'p'},
        {"no-gc-correction",     no_argument,       NULL, 0},
        {"num-samples",          required_argument, NULL, 'N'},
        {"burnin",               required_argument, NULL, 'B'},
        {"tss-cluster-distance", required_argument, NULL, 0},
        {0, 0, 0, 0}
    };

    Logger::level logger_level = Logger::INFO;
    unsigned int burnin = 100;
    unsigned int num_samples = 250;
    pos_t tss_cluster_dist = 30;
    bool run_gc_correction = true;
    constants::num_threads = boost::thread::hardware_concurrency();
    const char* fa_fn  = NULL;
    const char* output_filename = "isolator_output.h5";
    bool use_introns = false;
    bool use_exons = false;

    int opt;
    int optidx;
    while (true) {
        opt = getopt_long(argc, argv, "ho:vg:p:N:", long_options, &optidx);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_analyze_help(stdout);
                return 0;

            case 'o':
                output_filename = optarg;
                break;

            case 'v':
                logger_level = Logger::DEBUG;
                break;

            case 'p':
                constants::num_threads = std::max(1, atoi(optarg));
                break;

            case 'g':
                fa_fn = optarg;
                break;

            case 'B':
                burnin = strtoul(optarg, NULL, 10);
                break;

            case 'N':
                num_samples = strtoul(optarg, NULL, 10);
                break;

            case 0:
                if (strcmp(long_options[optidx].name, "no-gc-correction") == 0) {
                    run_gc_correction = false;
                }
                else if (strcmp(long_options[optidx].name, "tss-cluster-distance") == 0) {
                    tss_cluster_dist = strtod(optarg, NULL);
                }
                else if (strcmp(long_options[optidx].name, "introns") == 0) {
                    use_introns = true;
                }
                else if (strcmp(long_options[optidx].name, "exons") == 0) {
                    use_exons = true;
                }
                break;

            case '?':
                fprintf(stderr, "\n");
                print_analyze_help(stderr);
                return 1;

            default:
                abort();
        }
    }

    // no positional argumens
    if (optind + 1 >= argc) {
        print_analyze_usage(stdout);
        return 0;
    }

    print_logo();
    Logger::start();
    Logger::set_level(logger_level);

    // read transcripts
    const char* annotation_filename = argv[optind++];
    TranscriptSet ts;

    // TODO: parse a; bed file if the right option is set
    if (use_introns && use_exons) {
        Logger::abort("'--introns' and '--exons' arguments are mutually exclusize.");
    }
    else if (use_introns) {
        ts.read_bed(annotation_filename);
    }
    else if (use_exons) {
        Logger::abort("Reading exons from BED files is not yet supported.");
    }
    else {
        ts.read_gtf(annotation_filename, tss_cluster_dist);
    }

    hid_t output_file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (output_file_id < 0) {
        Logger::abort("Unable to open %s for writing.", output_filename);
    }

    // write metadata
    write_cassette_exon_data(output_file_id, ts);

    // initialize
    IsolatorMetadata metadata;
    Analyze analyze(burnin, num_samples, ts, fa_fn, run_gc_correction);
    int condition_num = 1;
    for (; optind < argc; ++optind) {
        char condition_name[100];
        snprintf(condition_name, sizeof(condition_name), "condition_%d", condition_num);

        const char* fn;
        for (fn = strtok(argv[optind], ","); fn; fn = strtok(NULL, ",")) {
            analyze.add_sample(condition_name, fn);
            metadata.sample_filenames.push_back(fn);
            metadata.sample_conditions.push_back(condition_name);
            Logger::debug("Adding '%s' to condition '%s'", fn, condition_name);
        }

        ++condition_num;
    }

    boost::timer::cpu_timer timer;
    timer.start();

    analyze.run(output_file_id);

    boost::timer::cpu_times elapsed_time = timer.elapsed();

    char elapsed_seconds_string[100];
    snprintf(elapsed_seconds_string, 100, "%0.1f", (double) elapsed_time.wall / 1e9);
    metadata.elapsed_seconds = elapsed_seconds_string;

    metadata.command_line = "isolator";
    for (int i = 0; i < argc; ++i) {
        metadata.command_line += " ";
        metadata.command_line += argv[i];
    }

    metadata.version = VERSION;

    time_t t = time(NULL);
    struct tm timestruct;
    localtime_r(&t, &timestruct);
    char time_string[200];
    strftime(time_string, 200, "%a, %d %b %Y %T %z", &timestruct);
    metadata.date = time_string;

    write_metadata(output_file_id, metadata);
    Logger::info("Finished. Have a nice day!");
    Logger::end();

    return EXIT_SUCCESS;
}


void print_usage(FILE* fout)
{
    fprintf(fout,
            "Usage: isolator <command> [<args>]\n"
            "Isolator version %s\n"
            "Instruction set: %s\n\n"
            "Where <command> is one of:\n"
            "    analyze           Quantify and test for differential expression\n"
            "                      and splicing, among other things.\n"
            "    summarize         Summarize a sampler run.\n"
            "    describe          Briefly describe contents of an output file.\n"
            "    help              Become enlightened.\n",
            VERSION, LINALG_INSTR_SET);
}


int isolator_help(int argc, char* argv[])
{
    if (argc <= 1) {
        print_usage(stdout);
        return EXIT_SUCCESS;
    }

    size_t cmdlen = 9 + strlen(argv[1]) + 1;
    char* cmd = new char [cmdlen];
    snprintf(cmd, cmdlen, "isolator-%s", argv[1]);

    execlp("man", "man", cmd, (char*) NULL);

    // if execlp returned, there was an error
    fprintf(stderr, "No suitable man page viewer found.\n");
    return EXIT_FAILURE;
}


int main(int argc, char* argv[])
{
    linalg_init();

    if (argc <= 1) {
        print_usage(stdout);
        return EXIT_SUCCESS;
    }

    ++argv;
    --argc;

    if (strcmp(argv[0], "-h") == 0 || strcmp(argv[0], "--help") == 0) {
        return isolator_help(argc, argv);
    }

    if (strcmp(argv[0], "summarize") == 0) {
        return isolator_summarize(argc, argv);
    }
    else if (strcmp(argv[0], "describe") == 0) {
        return isolator_describe(argc, argv);
    }
    else if (strcmp(argv[0], "analyze") == 0) {
        return isolator_analyze(argc, argv);
    }
    else if (strcmp(argv[0], "help") == 0) {
        return isolator_help(argc, argv);
    }

    fprintf(stderr, "Unknown command %s.\n\n", argv[0]);
    return EXIT_FAILURE;
}


