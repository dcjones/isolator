
#include "gtf_parse.h"
#include <string.h>
#include <ctype.h>
#include <math.h>


static const size_t init_str_size      = 128;
static const size_t gtf_file_buf_size  = 4096;

static void str_expand(str_t* s);
static void str_append(str_t* s, char c);
static void str_clear(str_t* s);

enum {
    STRAND_POS = 0,
    STRAND_NEG,
    STRAND_NA
};


static void* malloc_or_die(size_t n)
{
    void* ptr = malloc(n);
    if (ptr == NULL) {
        fprintf(stderr, "Failed to allocate %zu bytes. Out of memory!\n", n);
        exit(1);
    }

    return ptr;
}


static void* realloc_or_die(void* ptr, size_t n)
{
    ptr = realloc(ptr, n);
    if (ptr == NULL) {
        fprintf(stderr, "Failed to allocate %zu bytes. Out of memory!\n", n);
        exit(1);
    }

    return ptr;
}



static void fail(const char* msg, size_t line, size_t col)
{
    fprintf(stderr, "line: %zu, col: %zu\n", line, col);
    fprintf(stderr, "gtf parse error: %s\n", msg);
    exit(1);
}


str_t* str_alloc()
{
    str_t* s = malloc_or_die(sizeof(str_t));
    s->s = malloc_or_die(init_str_size);
    s->n = 0;
    s->size = init_str_size;
    str_clear(s);

    return s;
}


void str_free(str_t* s)
{
    free(s->s);
    free(s);
}


static void str_expand(str_t* s)
{
    s->size *= 2;
    s->s = realloc_or_die(s->s, s->size);
    memset(s->s + s->size / 2, '\0', s->size / 2);
}


static void str_append(str_t* s, char c)
{
    if (s->n + 1 >= s->size) str_expand(s);
    s->s[s->n++] = c;
}


static void str_clear(str_t* s)
{
    s->s[0] = '\0';
    s->n = 0;
}


void str_copy(str_t* dest, str_t* src)
{
    while (dest->size < src->n + 1) str_expand(dest);
    str_clear(dest);
    memcpy(dest->s, src->s, src->n + 1);
    dest->n = src->n;
}


gtf_row_t* gtf_row_alloc()
{
    gtf_row_t* row = malloc_or_die(sizeof(gtf_row_t));

    row->seqname = str_alloc();
    row->source  = str_alloc();
    row->feature = str_alloc();
    row->start      = -1;
    row->end        = -1;
    row->score      = 0.0;
    row->strand     = STRAND_NA;
    row->frame      = -1;
    row->attributes = str_map_create();

    return row;
}


static void gtf_row_clear(gtf_row_t* row)
{
    /* Attributes must be cleared because they may not be overwritten.
     * Everything else is, so we can leave it. */

    str_clear(row->seqname);
    str_clear(row->source);
    str_clear(row->feature);

    str_map_pair* u;
    size_t i;
    for (i = 0; i < row->attributes->n; ++i) {
        u = row->attributes->A[i];
        while (u) {
            str_clear((str_t*)u->value);
            u = u->next;
        }
    }
}


void gtf_row_free(gtf_row_t* row)
{
    str_free(row->seqname);
    str_free(row->source);
    str_free(row->feature);

    size_t i;
    str_map_pair* u;
    for (i = 0; i < row->attributes->n; ++i) {
        u = row->attributes->A[i];
        while (u) {
            str_free((str_t*)u->value);
            u = u->next;
        }
    }
    str_map_destroy(row->attributes);

    free(row);
}


typedef enum
{
    STATE_FIELD,      /* reading a field */
    STATE_KEY_SEEK,   /* seeking the next key */
    STATE_KEY,        /* reading a key */
    STATE_VALUE_SEEK, /* seeking the next value */
    STATE_VALUE,      /* reading a value */
    STATE_VALUE_END,
    STATE_EOF
} gtf_file_state;


gtf_file_t* gtf_file_alloc(FILE* file)
{
    gtf_file_t* f = malloc_or_die(sizeof(gtf_file_t));
    f->file   = file;
    f->k      = 0;
    f->state  = STATE_FIELD;
    f->buf    = malloc_or_die(gtf_file_buf_size);
    f->buf[0] = '\0';
    f->quote  = '\0';
    f->c      = f->buf;

    f->field1 = str_alloc();
    f->field2 = str_alloc();

    f->line = 1;
    f->col  = 1;

    return f;
}


void gtf_file_free(gtf_file_t* f)
{
    free(f->buf);
    str_free(f->field1);
    str_free(f->field2);
    free(f);
}


static void gtf_file_refill(gtf_file_t* f)
{
    char* s = fgets(f->buf, gtf_file_buf_size, f->file);

    if (s == NULL) {
        f->state = STATE_EOF;
        f->buf[0] = '\0';
    }

    f->c = f->buf;
}


bool gtf_next(gtf_file_t* f, gtf_row_t* r)
{
    gtf_row_clear(r);

    if (f->state == STATE_EOF) return false;

    int last_k = -1;
    str_t* field = NULL;
    char* endptr;
    str_t* attr;

    while (true) {
        /* figure out what string we are reading into */
        if (f->k != last_k) {
            last_k = f->k;
            switch (f->k) {
                case 0:
                    field = r->seqname;
                    break;

                case 1:
                    field = r->source;
                    break;

                case 2:
                    field = r->feature;
                    break;

                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                    field = f->field1;
                    break;

                default:
                    break;
            }
        }

        /* read more, if needed */
        if (*f->c == '\0') {
            gtf_file_refill(f);
            if (f->state == STATE_EOF) return false;
            continue;
        }

        if (*f->c == '\n') {
            /* skip blank lines */
            if (f->state == STATE_FIELD && f->k == 0 && field->n == 0) {
                f->c++;
                f->col = 1;
                f->line++;
            }
            else if (f->state != STATE_KEY_SEEK && f->state != STATE_VALUE_END &&
                     f->state != STATE_VALUE) {
                fail("Premature end of line.", f->line, f->col);
            }
        }

        /* handle the next character */
        switch (f->state) {
            case STATE_FIELD:
                if (*f->c == '\t') {
                    field->s[field->n] = '\0';

                    switch (f->k) {
                        case 3:
                            r->start = strtol(field->s, &endptr, 10);
                            if (*endptr != '\0') {
                                fail("Invalid start position.", f->line, f->col);
                            }
                            str_clear(field);
                            break;

                        case 4:
                            r->end = strtol(field->s, &endptr, 10);
                            if (*endptr != '\0') {
                                fail("Invalid end position.", f->line, f->col);
                            }
                            str_clear(field);
                            break;

                        case 5:
                            if (field->n == 1 && field->s[0] == '.') {
                                r->score = 0;
                            }
                            else {
                                r->score = strtod(field->s, &endptr);
                                if (*endptr != '\0') {
                                    fail("Invalid score.", f->line, f->col);
                                }
                            }
                            str_clear(field);
                            break;

                        case 6:
                            if (field->n > 0 && field->s[0] == '+') {
                                r->strand = STRAND_POS;
                            }
                            else if (field->n > 0 && field->s[0] == '-') {
                                r->strand = STRAND_NEG;
                            }
                            else {
                                r->strand = STRAND_NA;
                            }
                            str_clear(field);
                            break;

                        case 7:
                            if (field->n == 1 && field->s[0] == '.') {
                                r->frame = -1;
                            }
                            else {
                                r->frame = strtol(field->s, &endptr, 10);
                                if (*endptr != '\0') {
                                    fail("Invalid frame.", f->line, f->col);
                                }
                            }

                            str_clear(field);
                            break;

                        default:
                            break;
                    }

                    f->k++;
                    if (f->k == 8) f->state = STATE_KEY_SEEK;
                }
                else {
                    str_append(field, *f->c);
                }

                break;

            case STATE_KEY_SEEK:
                if (*f->c == '\n') {
                    f->c++;
                    goto gtf_next_finish;
                }

                if (isspace(*f->c)) break;

                f->state = STATE_KEY;

                if (*f->c == '"' || *f->c == '\'') {
                    f->quote = *f->c;
                    break;
                }
                else continue;

            case STATE_KEY:
                if ((f->quote != '\0' && *f->c == f->quote) ||
                    (f->quote == '\0' && isspace(*f->c))) {

                    f->field1->s[f->field1->n] = '\0';
                    f->state = STATE_VALUE_SEEK;
                    f->quote = '\0';
                }
                else {
                    str_append(f->field1, *f->c);
                }

                break;

            case STATE_VALUE_SEEK:
                if (isspace(*f->c)) break;

                f->state = STATE_VALUE;

                if (*f->c == '"' || *f->c == '\'') {
                    f->quote = *f->c;
                    break;
                }
                else continue;

            case STATE_VALUE:
                if ((*f->c == '\n') ||
                    (f->quote != '\0' && *f->c == f->quote) ||
                    (f->quote == '\0' && (*f->c == ';' || isspace(*f->c)))) {


                    if (*f->c == '\n' && f->quote != '\0') {
                        fail("Newline found before end quote.", f->line, f->col);
                    }

                    f->field2->s[f->field2->n] = '\0';
                    f->state = STATE_VALUE_END;

                    attr = (str_t*)str_map_get(r->attributes,
                                               f->field1->s,
                                               f->field1->n);

                    if (attr == NULL) {
                        str_map_set(r->attributes,
                                    f->field1->s,
                                    f->field1->n,
                                    str_alloc());

                        attr = (str_t*)str_map_get(r->attributes,
                                                   f->field1->s,
                                                   f->field1->n);
                    }

                    str_copy((str_t*)attr, f->field2);

                    str_clear(f->field1);
                    str_clear(f->field2);

                    if (f->quote != '\0') {
                        f->quote = '\0';
                        break;
                    }
                    else continue;
                }
                else {
                    str_append(f->field2, *f->c);
                    break;
                }


            case STATE_VALUE_END:
                if (*f->c != '\n' && isspace(*f->c)) break;
                else if (*f->c == ';') {
                    f->state = STATE_KEY_SEEK;
                    break;
                }
                else if (*f->c == '\n') {
                    f->c++;
                    goto gtf_next_finish;
                }
                else {
                    fail("Expected ';'.", f->line, f->col);
                }


            default:
                fputs("Inexplicable error in the gtf parser.\n", stderr);
                exit(1);
        }

        f->c++;
        f->col++;
    }

gtf_next_finish:

    f->state = STATE_FIELD;
    f->k     = 0;
    f->col   = 1;
    f->line++;
    return true;
}


