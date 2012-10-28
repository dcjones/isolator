
/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "murmurhash3.h"
#include "str_map.h"


static void* malloc_or_die(size_t n)
{
    void* ptr = malloc(n);
    if (ptr == NULL) {
        fprintf(stderr, "Failed to allocate %zu bytes. Out of memory!\n", n);
        exit(1);
    }

    return ptr;
}


static const size_t INITIAL_TABLE_SIZE = 16;
static const double MAX_LOAD = 0.77;


static void rehash(str_map* T, size_t new_n);
static void clear(str_map*);


str_map* str_map_create()
{
    str_map* T = malloc_or_die(sizeof(str_map));
    T->A = malloc_or_die(INITIAL_TABLE_SIZE * sizeof(str_map_pair*));
    memset(T->A, 0, INITIAL_TABLE_SIZE * sizeof(str_map_pair*));
    T->n = INITIAL_TABLE_SIZE;
    T->m = 0;
    T->max_m = T->n * MAX_LOAD;

    return T;
}


void str_map_destroy(str_map* T)
{
    if (T != NULL) {
        clear(T);
        free(T->A);
        free(T);
    }
}


void clear(str_map* T)
{
    str_map_pair* u;
    size_t i;
    for (i = 0; i < T->n; i++) {
        while (T->A[i]) {
            u = T->A[i]->next;
            free(T->A[i]->key);
            free(T->A[i]);
            T->A[i] = u;
        }
    }

    T->m = 0;
}


static void insert_without_copy(str_map* T, str_map_pair* V)
{
    uint32_t h = hash(V->key, V->keylen) % T->n;
    V->next = T->A[h];
    T->A[h] = V;
    T->m++;
}


static void rehash(str_map* T, size_t new_n)
{
    str_map U;
    U.n = new_n;
    U.m = 0;
    U.max_m = U.n * MAX_LOAD;
    U.A = malloc_or_die(U.n * sizeof(str_map_pair*));
    memset(U.A, 0, U.n * sizeof(str_map_pair*));

    str_map_pair *j, *k;
    size_t i;
    for (i = 0; i < T->n; i++) {
        j = T->A[i];
        while (j) {
            k = j->next;
            insert_without_copy(&U, j);
            j = k;
        }
        T->A[i] = NULL;
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = U.max_m;
}


void str_map_set(str_map* T, const char* key, size_t keylen, void* value)
{
    if (T->m >= T->max_m) rehash(T, T->n * 2);

    uint32_t h = hash(key, keylen) % T->n;

    str_map_pair* u = T->A[h];

    while (u) {
        if (u->keylen == keylen && memcmp(u->key, key, keylen) == 0) {
            u->value = value;
            return;
        }

        u = u->next;
    }

    u = malloc_or_die(sizeof(str_map_pair));
    u->key = malloc_or_die(keylen);
    memcpy(u->key, key, keylen);
    u->keylen = keylen;
    u->value  = value;

    u->next = T->A[h];
    T->A[h] = u;

    T->m++;
}


void* str_map_get(const str_map* T, const char* key, size_t keylen)
{
    uint32_t h = hash(key, keylen) % T->n;

    str_map_pair* u = T->A[h];

    while (u) {
        if (u->keylen == keylen && memcmp(u->key, key, keylen) == 0) {
            return u->value;
        }

        u = u->next;
    }

    return NULL;
}


void str_map_clear(str_map* T, void (*dealloc)(void*))
{
    size_t i;
    str_map_pair *u, *tmp;
    for (i = 0; i < T->n; ++i) {
        u = T->A[i];
        while (u != NULL) {
            dealloc(u->value);
            tmp = u;
            u = u->next;
            free(tmp);
        }

        T->A[i] = NULL;
    }
}

