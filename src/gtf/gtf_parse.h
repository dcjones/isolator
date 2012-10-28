/*
 * gtf_parse: parse the absolute shit out of GTF files.
 *
 * Daniel C. Jones <dcjones@cs.washington.edu>
 * 2011.04.29
 *
 */


#ifndef ISOLATOR_GTF_PARSE
#define ISOLATOR_GTF_PARSE

#if defined(__cplusplus)
extern "C" {
#endif

#include "str_map.h"
#include <stdio.h>
#include <stdbool.h>


/* a string that keeps track of its memory */

typedef struct
{
    char*  s;    /* null-terminated string */
    size_t n;    /* length of s */
    size_t size; /* bytes allocated for s */
} str_t;


str_t* str_alloc(void);
void   str_free(str_t* s);
void   str_copy(str_t* dest, str_t* src);


/* a gtf line */

typedef struct
{
    str_t*   seqname;
    str_t*   source;
    str_t*   feature;
    long     start;
    long     end;
    double   score;
    int      strand;
    int      frame;
    str_map* attributes;

} gtf_row_t;


gtf_row_t* gtf_row_alloc(void);
void gtf_row_free(gtf_row_t*);


/* a gtf file */

typedef struct
{
    FILE* file;
    int   k;     /* field number */
    int   state; /* parser state */
    char  quote; /* quotation mark that must be matched */
    char* buf;
    char* c;

    str_t* field1;
    str_t* field2;

    size_t line, col;

} gtf_file_t;


gtf_file_t* gtf_file_alloc(FILE* file);
void gtf_file_free(gtf_file_t*);
bool gtf_next(gtf_file_t*, gtf_row_t*);


#if defined(__cplusplus)
}
#endif

#endif



