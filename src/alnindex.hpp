
#ifndef ISOLATOR_ALNINDEX_HPP
#define ISOLATOR_ALNINDEX_HPP

#include <cstdlib>
#include <boost/thread.hpp>

#include "hat-trie/hat-trie.h"


/* Assign indexes to a set of string. (Read ids in this case.) */
class AlnIndex
{
    public:
        AlnIndex();
        ~AlnIndex();

        size_t size() const;
        void clear();

        long add(const char* key);

        /* Return -1 if the key is not present, otherwise return the key's
         * index. */
        long get(const char* key);

        size_t used_memory() const;

    private:
        hattrie_t* t;
        boost::mutex mut;
};

#endif

