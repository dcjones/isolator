/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#ifndef ISOLATOR_TABLE_H
#define ISOLATOR_TABLE_H

/**
 * \file
 *
 * A quick little hash table designed to be very good at one thing: hashing the
 * positions of every read in a BAM file.
 */


#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include "samtools/sam.h"


/* The table maps positions to counts. */
struct pos_table_val
{
    int32_t  pos;
    uint32_t count;
};


/* Each strand of each sequence is stored in a seperate hash table, since the
 * BAM file already assigns an integer index to sequences */
struct pos_subtable
{
    struct pos_table_val* A; /* table proper */
    size_t n;               /* table size (as an index into a prime table) */
    size_t m;               /* hashed items */
    size_t max_m;           /* max hashed items before rehash */
};


struct pos_table
{
    /* an array indexd by strand -> sequence id */
    struct pos_subtable* ts[2];
    size_t m; /* number of unique positions */
    size_t n; /* number of sequences */
    char** seq_names;
};



/* initialize, where n is the number of sequences */
void pos_table_create(struct pos_table* T, size_t n);
void pos_table_copy(struct pos_table* T, const struct pos_table* U);
void pos_table_destroy(struct pos_table* T);

void pos_table_inc(struct pos_table*, bam1_t* read);
void pos_table_inc_pos(struct pos_table*, int32_t tid, int32_t pos, uint32_t strand);

uint32_t pos_table_count(struct pos_table*, bam1_t* read);
uint32_t pos_table_count_pos(struct pos_table*, int32_t tid, int32_t pos, uint32_t strand);





/* A simple structure storing the read count at each nonzero position across the
 * genome, sorted by position to allow for binary search. Using this we can
 * get the read count across any interval very quickly. */
struct read_counts
{
    struct pos_table_val** xss[2]; /* read count arrays indexed by strand -> tid */
    size_t* mss[2]; /* lengths of xss arrays */
    size_t m;       /* total number of unique positions */
    size_t n;       /* number of sequences */
    char**  seq_names;
};


void read_counts_create( struct read_counts* C, const struct pos_table* T );
void read_counts_copy( struct read_counts* C, const struct read_counts* B );
void read_counts_destroy( struct read_counts* C );


void read_counts_count( const struct read_counts* C,
                        int32_t tid, int32_t start, int32_t end, uint32_t strand,
                        unsigned int* ys );

unsigned int read_counts_total( const struct read_counts* C,
                                int32_t tid, int32_t start, int32_t end, uint32_t strand );


void read_count_occurances( const struct read_counts* C,
                            int32_t tid,  int32_t start, int32_t end, uint32_t strand,
                            uint64_t* ks, size_t max_k );



/* dump to a flat array (used by sequencing bias) */
struct read_pos
{
    int32_t  tid;
    uint32_t strand;
    int32_t  pos;
    uint32_t count;
};

void pos_table_dump( struct pos_table* T, struct read_pos** A_, size_t* N_, size_t limit );

#ifdef __cplusplus
}
#endif

#endif




