
#include <cstdlib>
#include <getopt.h>
#include <unistd.h>

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
            "-h, --help                print this help message\n"
            "-o, --out=FILE            output results to the given file (default: fpkm.tab)\n"
            "-g, --genomic-seq=FILE    correct for sequence bias, given the a the sequence\n"
            "                          against which the reads are aligned, in FAST format.\n"
            "-p, --threads=N           number of threads to use.\n"
            "-T, --trans-ids=FILE      a file listing transcript id's of transcripts that\n"
            "                          will be quantified (by default, every transcript).\n"
            "-G, --gene-ids=FILE       a file listing gene id's of transcripts that\n"
            "                          will be quantified (by default, every transcript).\n\n"
            "See 'isolator help quantify' for more.\n");
}


int quantify(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"help",        no_argument,       NULL, 'h'},
        {"out",         required_argument, NULL, 'o'},
        {"threads",     required_argument, NULL, 'p'},
        {"trans-ids",   required_argument, NULL, 'T'},
        {"gene-ids",    required_argument, NULL, 'G'},
        {"genomic-seq", required_argument, NULL, 'g'},
        {0, 0, 0, 0}
    };

    const char* fa_fn  = NULL;
    const char* out_fn = "isolator.db";
    constants::num_threads = boost::thread::hardware_concurrency();

    int opt;
    int opt_idx;

    while (true) {
        opt = getopt_long(argc, argv, "ho:p:g:", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'h':
                print_quantify_help(stdout);
                return 0;

            case 'o':
                out_fn = optarg;
                break;

            case 'p':
                constants::num_threads = std::max(1, atoi(optarg));
                break;

            case 'g':
                fa_fn = optarg;
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
    else if (optind + 2 < argc) {
        fprintf(stderr, "Too many arguments.\n\n");
        print_quantify_usage(stderr);
        return 1;
    }

    print_logo();
    Logger::start();

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
    FragmentModel fm;
    fm.estimate(ts, bam_fn, fa_fn);

    /* Initialize the sampler. */
    Sampler sampler(bam_fn, fa_fn, ts, fm);
    /* TODO: run the sampler */

    Logger::info("Finished. Have a nice day!");
    Logger::end();

    return EXIT_SUCCESS;
}


void print_usage(FILE* fout)
{
    fprintf(fout,
            "Usage: isolator <command> [<args>]\n"
            "Isolator version %s\n\n"
            "Where <command> is one of:\n"
            "    quantify          Quantify transcript abundance.\n"
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
    else if (strcmp(argv[0], "test") == 0) {
        return EXIT_FAILURE; /* TODO */
    }
    else if (strcmp(argv[0], "help") == 0) {
        return isolator_help(argc, argv);
    }

    fprintf(stderr, "Unknown command %s.\n\n", argv[0]);
    return EXIT_FAILURE;
}


