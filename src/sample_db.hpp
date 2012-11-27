
#ifndef ISOLATOR_SAMPLE_DB
#define ISOLATOR_SAMPLE_DB

#include <boost/iterator/iterator_facade.hpp>

#include "common.hpp"
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

        void begin_transaction();
        void commit_transaction();

        void insert_meta(const char* key, const char* val);
        void insert_param(const char* key, double param);

        void insert_sampler_result(TranscriptID transcript_id,
                                   GeneID gene_id,
                                   float effective_length,
                                   float map_estimate,
                                   float* samples,
                                   size_t num_samples);

    private:
        void exec(const char* stmt);
        sqlite3_stmt* prep(const char* stmt);

        sqlite3* db;

        /* prepared statements */
        sqlite3_stmt* meta_ins_stmt;
        sqlite3_stmt* param_ins_stmt;
        sqlite3_stmt* sample_ins_stmt;
        sqlite3_stmt* id_ins_stmt;

        friend class SampleDBIterator;
};


struct SampleDBEntry
{
    TranscriptID transcript_id;
    GeneID gene_id;
    float effective_length;
    float map_estimate;
    std::vector<float> samples;
};


/* Iterate through the sampler output for each transcript. */
class SampleDBIterator :
    public boost::iterator_facade<SampleDBIterator,
                                  const SampleDBEntry,
                                  boost::forward_traversal_tag>
{
    public:
        SampleDBIterator();
        SampleDBIterator(SampleDB&);
        ~SampleDBIterator();

    private:
        friend class boost::iterator_core_access;

        void increment();
        bool equal(const SampleDBIterator& other) const;
        const SampleDBEntry& dereference() const;

        SampleDB* owner;
        sqlite3_stmt* s;
        unsigned int num_samples;

        SampleDBEntry entry;
};


#endif

