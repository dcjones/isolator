
SUBDIRS = gtf hat-trie samtools seqbias sqlite3

bin_PROGRAMS = isolator

isolator_SOURCES = isolator.cpp \
                   alloc.hpp \
                   common.hpp common.cpp \
                   constants.hpp constants.cpp \
                   emp_dist.hpp emp_dist.cpp \
                   fragment_model.hpp fragment_model.cpp \
                   intervals.hpp intervals.cpp \
                   linalg.hpp linalg.cpp \
                   logger.hpp logger.cpp \
                   pos_table.hpp pos_table.cpp \
                   queue.hpp \
                   read_set.hpp read_set.cpp \
                   sampler.hpp sampler.cpp \
                   sample_db.hpp sample_db.cpp \
                   sparse_mat.hpp sparse_mat.cpp \
                   transcripts.hpp transcripts.cpp \
                   trie.hpp

isolator_LDADD = gtf/libgtf.la \
                 hat-trie/libhat-trie.la \
                 samtools/libbam.la \
                 seqbias/libseqbias.la \
                 sqlite3/libsqlite3.la \
                 $(BOOST_THREAD_LIBS) \
                 $(BOOST_SYSTEM_LIBS) \
                 $(BOOST_TIMER_LIBS)

noinst_PROGRAMS = hattrie-bench

hattrie_bench_SOURCES = hattrie-bench.c

hattrie_bench_LDADD = hat-trie/libhat-trie.la \
                      samtools/libbam.la \
                      -lz

