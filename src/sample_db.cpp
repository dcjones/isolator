
#include <zlib.h>

#include "sample_db.hpp"
#include "logger.hpp"

SampleDB::SampleDB(const char* fn, bool writeable)
{
    int flags = writeable ? SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE :
                            SQLITE_OPEN_READONLY;

    int result = sqlite3_open_v2(fn, &db, flags, NULL);
    if (result != SQLITE_OK) {
        Logger::abort("Could not open sqlite3 database: '%s'.",
                      sqlite3_errmsg(db));
    }

    if (flags & SQLITE_OPEN_CREATE) {
        /* create schemea */
        exec("drop table if exists meta");
        exec("create table meta (key text, value text)");

        exec("drop table if exists ids");
        exec("create table ids (transcript_id text primary key, gene_id text)");

        exec("drop table if exists parameters");
        exec("create table parameters (name text, value real)");

        exec("drop table if exists samples");
        exec("create table samples (transcript_id text primary key, effective_length real, samples blob)");
    }

    /* prepared statements */
    meta_ins_stmt   = prep("insert into meta values (?, ?)");
    param_ins_stmt  = prep("insert into parameters values (?, ?)");
    sample_ins_stmt = prep("insert into samples values (?, ?, ?, ?)");
    id_ins_stmt     = prep("insert into ids values (?, ?)");
}


SampleDB::~SampleDB()
{
    sqlite3_finalize(meta_ins_stmt);
    sqlite3_finalize(param_ins_stmt);
    sqlite3_finalize(sample_ins_stmt);
    sqlite3_finalize(id_ins_stmt);
    sqlite3_close(db);
}


void SampleDB::begin_transaction()
{
    exec("begin transaction");
}


void SampleDB::commit_transaction()
{
    exec("commit transaction");
}


void SampleDB::insert_meta(const char* key, const char* val)
{
    sqlite3_bind_text(meta_ins_stmt, 1, key, -1, SQLITE_STATIC);
    sqlite3_bind_text(meta_ins_stmt, 2, val, -1, SQLITE_STATIC);
    int result = sqlite3_step(meta_ins_stmt);
    if (result != SQLITE_DONE) {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(db));
    }
    sqlite3_reset(meta_ins_stmt);
}


void SampleDB::insert_param(const char* key, double val)
{
    sqlite3_bind_text(param_ins_stmt, 1, key, -1, SQLITE_STATIC);
    sqlite3_bind_double(param_ins_stmt, 2, val);
    int result = sqlite3_step(param_ins_stmt);
    if (result != SQLITE_DONE) {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(db));
    }
    sqlite3_reset(param_ins_stmt);
}


void SampleDB::insert_sampler_result(TranscriptID transcript_id,
                                     GeneID gene_id,
                                     float effective_length,
                                     float* samples,
                                     size_t num_samples)
{

    sqlite3_bind_text(sample_ins_stmt, 1, transcript_id.get().c_str(),
                      -1, SQLITE_STATIC);
    sqlite3_bind_double(sample_ins_stmt, 2, effective_length);
    Bytef* zbuf = NULL;

    size_t N  = num_samples * sizeof(float);
    uLong bound = compressBound(N);
    zbuf = new Bytef [bound];
    uLong readsize = bound;
    int ret = compress2(zbuf, &readsize,
                        reinterpret_cast<const Bytef*>(samples),
                        N, Z_BEST_COMPRESSION);
    if (ret != Z_OK) {
        Logger::abort("Error compressing results: %d.", ret);
    }

    sqlite3_bind_blob(sample_ins_stmt, 3,
                     reinterpret_cast<const void*>(zbuf),
                     (int) readsize,
                     SQLITE_STATIC);
    int result = sqlite3_step(sample_ins_stmt);

    if (result != SQLITE_DONE) {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(db));
    }

    delete [] zbuf;

    sqlite3_reset(sample_ins_stmt);

    /* insert gene_id / transcript_id */
    sqlite3_bind_text(id_ins_stmt, 1, transcript_id.get().c_str(),
                      -1, SQLITE_STATIC);
    sqlite3_bind_text(id_ins_stmt, 2, gene_id.get().c_str(),
                      -1, SQLITE_STATIC);

    result = sqlite3_step(id_ins_stmt);

    if (result != SQLITE_DONE) {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(db));
    }

    sqlite3_reset(id_ins_stmt);
}


void SampleDB::exec(const char* stmt)
{
    int result;
    char* errmsg;
    result = sqlite3_exec(db, stmt, NULL, NULL, &errmsg);

    if (result != SQLITE_OK) {
        Logger::abort("Sqlite3 error: '%s'.", errmsg);
    }
}


sqlite3_stmt* SampleDB::prep(const char* stmt)
{
    sqlite3_stmt* s;
    int result = sqlite3_prepare_v2(db, stmt, -1, &s, NULL);

    if (result != SQLITE_OK) {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(db));
    }

    return s;
}


SampleDBIterator::SampleDBIterator()
    : owner(NULL)
    , s(NULL)
{
}


SampleDBIterator::SampleDBIterator(SampleDB& owner)
    : owner(&owner)
{
    int result = sqlite3_prepare_v2(owner.db,
            "select value from parameters where name = \"num_samples\"", -1, &s, NULL);
    if (result != SQLITE_OK) {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(owner.db));
    }

    result = sqlite3_step(s);
    if (result != SQLITE_ROW) {
        Logger::abort("Malformed sampler results.");
    }
    num_samples = (unsigned int) sqlite3_column_double(s, 0);
    sqlite3_finalize(s);
    entry.samples.resize(num_samples);

    result = sqlite3_prepare_v2(owner.db,
            "select transcript_id, gene_id, effective_length, samples "
            "from samples join ids using (transcript_id)", -1, &s, NULL);

    if (result != SQLITE_OK) {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(owner.db));
    }
    increment();
}


SampleDBIterator::~SampleDBIterator()
{
    if (s) sqlite3_finalize(s);
}


void SampleDBIterator::increment()
{
    if (owner == NULL || s == NULL) return;

    int result = sqlite3_step(s);
    if (result == SQLITE_ROW) {
        entry.transcript_id    = (const char*) sqlite3_column_text(s, 0);
        entry.gene_id          = (const char*) sqlite3_column_text(s, 1);
        entry.effective_length = (float) sqlite3_column_double(s, 2);

        const void* blob = sqlite3_column_blob(s, 3);
        size_t blobsize = sqlite3_column_bytes(s, 3);

        uLongf N = num_samples * sizeof(float);
        int ret = uncompress(reinterpret_cast<Bytef*>(&entry.samples.at(0)), &N,
                             reinterpret_cast<const Bytef*>(blob), blobsize);
        if (ret != Z_OK) {
            Logger::abort("Error decompressing sampler output %d.", ret);
        }

    }
    else if (result == SQLITE_DONE) s = NULL;
    else {
        Logger::abort("Sqlite3 error: '%s'.", sqlite3_errmsg(owner->db));
    }
}


bool SampleDBIterator::equal(const SampleDBIterator& other) const
{
    return s == other.s;
}


const SampleDBEntry& SampleDBIterator::dereference() const
{
    return entry;
}



