
#ifndef ISOLATOR_TRIE_HPP
#define ISOLATOR_TRIE_HPP

#include <cstring>
#include <boost/iterator/iterator_facade.hpp>

#include "hat-trie/hat-trie.h"


/* These are wrappers for the C hat-trie implementation. It is a memory
 * efficient string table, particularly when the string stored tend to share
 * prefixes and order is not impartant. */


/* A set of strings. */
class TrieSet
{
    public:
        TrieSet()
        {
            t = hattrie_create();
        }

        ~TrieSet()
        {
            hattrie_free(t);
        }

        void add(const char* key)
        {
            *hattrie_get(t, key, strlen(key) + 1) = 1;
        }

        bool has(const char* key)
        {
            return hattrie_tryget(t, key, strlen(key) + 1) != NULL;
        }

    private:
        hattrie_t* t;
        friend class TrieSetIterator;
};


/* Unordered iteration over items in a TrieSet. */
class TrieSetIterator :
    public boost::iterator_facade<TrieSetIterator,
                                  const char*,
                                  boost::forward_traversal_tag>
{
    public:
        TrieSetIterator()
            : it(NULL)
        {}

        TrieSetIterator(const TrieSet& s)
        {
            it = hattrie_iter_begin(s.t, false);
        }

        ~TrieSetIterator()
        {
            hattrie_iter_free(it);
        }

    private:
        friend class boost::iterator_core_access;

        void increment()
        {
            hattrie_iter_next(it);
        }

        bool equal(TrieSetIterator const& other) const
        {
            return ((it == NULL || hattrie_iter_finished(it)) &&
                    (other.it == NULL || hattrie_iter_finished(other.it))) ||
                   hattrie_iter_equal(it, other.it);
        }

        const char* dereference() const
        {
            size_t len;
            return hattrie_iter_key(it, &len);
        }

        hattrie_iter_t* it;
};


/* A map of string to an arbitrary object. */
template <typename T>
class TrieMap
{
    public:
        TrieMap()
        {
            t = hattrie_create();
        }

        ~TrieMap()
        {
            hattrie_free(t);
        }

        bool has(const char* key)
        {
            hattrie_tryget(t, key, strlen(key)) != NULL;
        }

        T& operator [] (const char* key)
        {
            T** slot = reinterpret_cast<T**>(
                            hattrie_get(t, key, strlen(key) + 1));
            if (*slot == NULL) {
                *slot = new T();
            }

            return **slot;
        }

    private:
        hattrie_t* t;

        template <typename U>
        friend class TrieMapIterator;
};


template <typename T>
class TrieMapIterator :
    public boost::iterator_facade<TrieMapIterator<T>,
                                  const std::pair<const char*, T&>,
                                  boost::forward_traversal_tag>
{
    public:

        TrieMapIterator()
            : it(NULL)
        {
        }

        TrieMapIterator(TrieMap<T>& s)
        {
            it = hattrie_iter_begin(s.t, false);
        }

        ~TrieMapIterator()
        {
            hattrie_iter_free(it);
        }


    private:
        friend class boost::iterator_core_access;

        void increment()
        {
            hattrie_iter_next(it);
        }

        bool equal(TrieMapIterator<T> const& other) const
        {
            return ((it == NULL || hattrie_iter_finished(it)) &&
                    (other.it == NULL || hattrie_iter_finished(other.it))) ||
                   hattrie_iter_equal(it, other.it);
        }

        const std::pair<const char*, T&> dereference() const
        {
            size_t len;
            const char* key = hattrie_iter_key(it, &len);
            T** slot = reinterpret_cast<T**>(hattrie_iter_val(it));
            return std::pair<const char*, T&>(key, **slot);
        }

        hattrie_iter_t* it;
};


#endif

