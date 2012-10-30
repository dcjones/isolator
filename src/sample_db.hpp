
#ifndef ISOLATOR_SAMPLE_DB
#define ISOLATOR_SAMPLE_DB

#include "sqlite3/sqlite3.h"

/* A wrapper for a sqlite3 database containing output from a sampler (i.e.,
 * 'isolator quantify') run. */
class SampleDB
{
    public:
        /* Open a SampleDB sqlite3 database.
         *
         * Args:
         *   fn: Filename.
         *   writeable: If true open for writing and create if nonexistant.
         */
        SampleDB(const char* fn, bool writable);
        ~SampleDB();

    private:
        void exec(const char* stmt);
        sqlite3_stmt* prep(const char* stmt);

        sqlite3* db;

        /* prepared statements */
        sqlite3_stmt* meta_ins_stmt;
        sqlite3_stmt* param_ins_stmt;
};

#endif

