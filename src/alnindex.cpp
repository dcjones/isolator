
#include "alnindex.hpp"
#include <cstring>


AlnIndex::AlnIndex()
{
    t = hattrie_create();
}


AlnIndex::~AlnIndex()
{
    hattrie_free(t);
}


size_t AlnIndex::size() const
{
    return hattrie_size(t);
}


void AlnIndex::clear()
{
    hattrie_clear(t);
}


long AlnIndex::add(const char* key)
{
    value_t* val = hattrie_get(t, key, strlen(key));
    if (*val == 0) {
        *val = hattrie_size(t);
    }
    return *val;
}


long AlnIndex::get(const char* key)
{
    value_t* val = hattrie_tryget(t, key, strlen(key));
    if (val == NULL) return -1;
    else {
        return *val - 1;;
    }
}


size_t AlnIndex::used_memory() const
{
    return hattrie_sizeof(t);
}


