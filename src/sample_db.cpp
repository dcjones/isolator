
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
        exec("create table samples (transcript_id text primary key, map_estimate real, samples blob)");
    }

    /* prepared statements */
    meta_ins_stmt   = prep("insert into meta values (?, ?)");
    param_ins_stmt  = prep("insert into parameters values (?, ?)");
}


SampleDB::~SampleDB()
{
    sqlite3_finalize(meta_ins_stmt);
    sqlite3_finalize(param_ins_stmt);
    sqlite3_close(db);
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


