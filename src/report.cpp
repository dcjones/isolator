
#include <hdf5.h>
#include <hdf5_hl.h>

#include "logger.hpp"
#include "report.hpp"
#include "summarize.hpp"


void generate_report(const char* input_filename, FILE* out)
{
    Summarize summarize(input_filename);
    IsolatorMetadata metadata;
    summarize.read_metadata(metadata);

    fprintf(out,
        "<!DOCTYPE html>\n"
        "<html>\n"
        "<head>\n"
        "  <title>Isolator Report</title>\n"
        "</head>\n"
        "<body>\n"
        "  <h1>Isolator Report</h1>\n");

    // Output basic information about the run:

    fprintf(out, "  <p>Isolator version: %s</p>\n", metadata.version.c_str());
    fprintf(out, "  <p>Run started at %s, finished in %s seconds.</p>\n",
            metadata.date.c_str(), metadata.elapsed_seconds.c_str());
    fprintf(out, "  <p>Command line: <pre>%s</pre></p>\n",
            metadata.command_line.c_str());

    fprintf(out,
        "  <h3>Samples</h3>\n"
        "  <table>\n"
        "    <tr>\n"
        "      <td>Sample Number</td> <td>Condition</td> <td>Sample Filename</td>\n"
        "    </tr>\n");

    for (size_t i = 0; i < metadata.sample_filenames.size(); ++i) {
        fprintf(out,
            "    <tr>\n"
            "      <td>%d</td> <td>%s</td> <td>%s</td>\n"
            "    </tr>\n",
            (int) i + 1,
            metadata.sample_conditions[i].c_str(),
            metadata.sample_filenames[i].c_str());
    }

    fprintf(out, "  </table>\n");

    fprintf(out,
        "  <script>\n"
        "    var isolator_data = {};\n");

       // dump gene_ids
    fprintf(out,
        "    isolator_data.gene_ids = [\n");


           //     ]

           // "  </script>\n");


       /* Ok, then. What is the raw data are we going to include?

       Tgroup expression:

       gene_id, gene_name, transcript_id, lower, mid, upper

       for every sample.

       */




    /* I'm kind of thinking about writing the raw data to json, embedding that
     * in the document then having some javascript that's responsible for generating
     * tables and plots.
     */

    // Summarize differential transcription
    fprintf(out, "\n\n  <h2>Differential Transcription</h2>\n");
    // TODO

    // Summarize differential splicing
    fprintf(out, "\n\n  <h2>Differential Splicing</h2>\n");
    // TODO

    fprintf(out,
        "</body>\n"
        "</html>\n");
}