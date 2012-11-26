
#include <cstdlib>
#include <getopt.h>
#include <unistd.h>
#include <gsl/gsl_statistics_float.h>

#include "config.h"
#include "constants.hpp"
#include "fragment_model.hpp"
#include "logger.hpp"
#include "sample_db.hpp"
#include "sampler.hpp"
#include "transcripts.hpp"

const char* isolator_logo =
"     _           _       _\n"
"    (_)___  ___ | | __ _| |_ ___  _ __\n"
"    | / __|/ _ \\| |/ _` | __/ _ \\| '__|\n"
"    | \\__ \\ (_) | | (_| | || (_) | |\n"
"    |_|___/\\___/|_|\\__,_|\\__\\___/|_|\n";


static void print_logo()
{
    printf("%s\n     Version: %s\n\n", isolator_logo, VERSION);
}


void print_quantify_usage(FILE* fout)
{
    fprintf(fout, "Usage: isolator quantify [options] genes.gtf reads.bam\n");
}


void print_quantify_help(FILE* fout)
{
    print_quantify_usage(fout);
    fprintf(fout,
            "\nOptions:\n"
            "-h, --help                Print this help message\n"
            "-v, --verbose             Print a bunch of information useful mainly for debugging.\n"
            "-o, --out=FILE            Output results to the given file (default: fpkm.tab)\n"
            "-g, --genomic-seq=FILE    Correct for sequence bias, given the a the sequence\n"
            "                          against which the reads are aligned, in FAST format.\n"
            "-p, --threads=N           number of threads to use.\n"
            "-T, --trans-ids=FILE      A file listing transcript id's of transcripts that\n"
            "                          will be quantified (by default, every transcript).\n"
            "-G, --gene-ids=FILE       A file listing gene id's of transcripts that\n"
            "                          will be quantified (by default, every transcript).\n"
            "-N, --num-samples         Generate this number of samples (by default: 10000).\n\n"
            "See 'isolator help quantify' for more.\n");
}


int quantify(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"help",        no_argument,       NULL, 'h'},
        {"verbose",     no_argument,       NULL, 'v'},
        {"out",         required_argument, NULL, 'o'},
        {"threads",     required_argument, NULL, 'p'},
        {"trans-ids",   required_argument, NULL, 'T'},
        {"gene-ids",    required_argument, NULL, 'G'},
        {"genomic-seq", required_argument, NULL, 'g'},
        {"num-samples", required_argument, NULL, 'N'},
        {0, 0, 0, 0}
    };

    const char* fa_fn  = NULL;
    const char* out_fn = "isolator.db";
    constants::num_threads = boost::thread::hardware_concurrency();
    unsigned int num_samples = 1000;
    Logger::level logger_level = Logger::INFO;

    int opt;
    int opt_idx;

    while (true) {
        opt = getopt_long(argc, argv, "hvo:p:g:N:", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_quantify_help(stdout);
                return 0;

            case 'v':
                logger_level = Logger::DEBUG;
                break;

            case 'o':
                out_fn = optarg;
                break;

            case 'p':
                constants::num_threads = std::max(1, atoi(optarg));
                break;

            case 'g':
                fa_fn = optarg;
                break;

            case 'N':
                num_samples = strtoul(optarg, NULL, 10);
                break;

            case '?':
                fprintf(stderr, "\n");
                print_quantify_help(stderr);
                return 1;

            default:
                abort();
        }
    }

    /* no positional argumens */
    if (optind == argc) {
        print_quantify_usage(stdout);
        return 0;
    }

    /* too few positional arguments */
    else if (optind + 1 == argc) {
        fprintf(stderr, "Too few arguments.\n\n");
        print_quantify_usage(stderr);
        return 1;
    }

    /* too many */
    else if (optind + 2 > argc) {
        fprintf(stderr, "Too many arguments.\n\n");
        print_quantify_usage(stderr);
        return 1;
    }

    print_logo();
    Logger::start();
    Logger::set_level(logger_level);

    const char* gtf_fn = argv[optind];
    const char* bam_fn = argv[optind + 1];

    /* Read transcripts. */
    TranscriptSet ts;
    FILE* gtf_f = fopen(gtf_fn, "rb");
    if (gtf_f == NULL) {
        Logger::abort("Can't open file %s for reading.", gtf_fn);
    }
    ts.read_gtf(gtf_f);
    fclose(gtf_f);

    /* Prepare output database. */
    SampleDB db(out_fn, true);

    /* Initialize the fragment model. */
    FragmentModel* fm = new FragmentModel();
    fm->estimate(ts, bam_fn, fa_fn);

    /* Initialize the sampler. */
    Sampler sampler(bam_fn, fa_fn, ts, *fm);
    delete fm; /* free a little memory */
    sampler.run(num_samples, db);

    Logger::info("Finished. Have a nice day!");
    Logger::end();

    return EXIT_SUCCESS;
}


void print_summarize_usage(FILE* fout)
{
    fprintf(fout, "Usage: isolator summarize [options] quantify_output.db\n");
}


void print_summarize_help(FILE* fout)
{
    print_summarize_usage(fout);
    fprintf(fout,
            "\nOptions:\n"
            "-h, --help                Print this help message\n"
            "-v, --verbose             Print a bunch of information useful mainly for debugging.\n"
            "-o, --out=FILE            Output summary to the given file. (default: standard out).\n"
            "\n"
            "See 'isolator help summarize' for more.\n");
}


int summarize(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"help",        no_argument,       NULL, 'h'},
        {"verbose",     no_argument,       NULL, 'v'},
        {"out",         required_argument, NULL, 'o'},
        {0, 0, 0, 0}
    };

    const char* out_fn = NULL;
    Logger::level logger_level = Logger::INFO;

    int opt;
    int opt_idx;

    while (true) {
        opt = getopt_long(argc, argv, "hvo:", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_quantify_help(stdout);
                return 0;

            case 'v':
                logger_level = Logger::DEBUG;
                break;

            case 'o':
                out_fn = optarg;
                break;

            case '?':
                fprintf(stderr, "\n");
                print_quantify_help(stderr);
                return 1;

            default:
                abort();
        }
    }

    /* no positional argumens */
    if (optind == argc) {
        print_quantify_usage(stdout);
        return 0;
    }

    /* too many */
    else if (optind + 1 > argc) {
        fprintf(stderr, "Too many arguments.\n\n");
        print_quantify_usage(stderr);
        return 1;
    }

    FILE* out_f = stdout;
    if (out_fn) {
        out_f = fopen(out_fn, "w");
        if (out_f == NULL) {
            fprintf(stderr, "Can't open file %s for writing.\n", out_fn);
        }
    }

    const char* in_fn = argv[optind];
    SampleDB db(in_fn, false);

    fprintf(out_f,
            "transcript_id\tgene_id\tmap_estimate\teffective_length\t"
            "posterior_mean\tposterior_sd\tposterior_median\tposterior_mad\t"
            "lower_95_cred\tupper_95_cred\n");

    std::vector<float> samples;
    for (SampleDBIterator i(db); i != SampleDBIterator(); ++i) {
        samples = i->samples;
        std::sort(samples.begin(), samples.end());
        float posterior_mean = gsl_stats_float_mean(
                &samples.at(0), 1, samples.size());
        float posterior_sd = gsl_stats_float_sd_m(
                &samples.at(0), 1, samples.size(), posterior_mean);
        float posterior_median = gsl_stats_float_median_from_sorted_data(
                &samples.at(0), 1, samples.size());

        float lower_95_cred = gsl_stats_float_quantile_from_sorted_data(
                &samples.at(0), 1, samples.size(), 0.025);
        float upper_95_cred = gsl_stats_float_quantile_from_sorted_data(
                &samples.at(0), 1, samples.size(), 0.975);

        for (std::vector<float>::iterator j = samples.begin();
                j != samples.end(); ++j) {
            *j = fabsf(*j - posterior_median);
        }
        std::sort(samples.begin(), samples.end());
        float posterior_mad = gsl_stats_float_median_from_sorted_data(
                &samples.at(0), 1, samples.size());

        fprintf(out_f, "%s\t%s\t%e\t%f\t%e\t%e\t%e\t%e\t%e\t%e\n",
                i->transcript_id.get().c_str(),
                i->gene_id.get().c_str(),
                i->map_estimate,
                i->effective_length,
                posterior_mean,
                posterior_sd,
                posterior_median,
                posterior_mad,
                lower_95_cred,
                upper_95_cred);
    }

    if (out_fn) fclose(out_f);
    return EXIT_SUCCESS;
}


void print_usage(FILE* fout)
{
    fprintf(fout,
            "Usage: isolator <command> [<args>]\n"
            "Isolator version %s\n\n"
            "Where <command> is one of:\n"
            "    quantify          Quantify transcript abundance.\n"
            "    summarize         Summarize a sampler run.\n"
            "    report            Generate useful output from an quantify run.\n"
            "    test              Test for differential expression and\n"
            "                      isoform switching.\n"
            "    help              Become enlightened.\n",
            VERSION);
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
    if (argc <= 1) {
        print_usage(stdout);
        return EXIT_SUCCESS;
    }

    ++argv;
    --argc;

    if (strcmp(argv[0], "-h") == 0 || strcmp(argv[0], "--help") == 0) {
        return isolator_help(argc, argv);
    }

    if (strcmp(argv[0], "quantify") == 0) {
        return quantify(argc, argv);
    }
    else if (strcmp(argv[0], "summarize") == 0) {
        return summarize(argc, argv);
    }
    else if (strcmp(argv[0], "test") == 0) {
        return EXIT_FAILURE; /* TODO */
    }
    else if (strcmp(argv[0], "help") == 0) {
        return isolator_help(argc, argv);
    }

    fprintf(stderr, "Unknown command %s.\n\n", argv[0]);
    return EXIT_FAILURE;
}


