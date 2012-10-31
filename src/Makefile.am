
SUBDIRS = gtf hat-trie samtools sqlite3

bin_PROGRAMS = isolator

isolator_SOURCES = isolator.cpp \
                   common.hpp common.cpp \
                   constants.hpp constants.cpp \
                   fragment_model.hpp fragment_model.cpp \
                   intervals.hpp intervals.cpp \
                   logger.hpp logger.cpp \
                   queue.hpp \
                   read_set.hpp read_set.cpp \
                   sam_scan.hpp sam_scan.cpp \
                   sample_db.hpp sample_db.cpp \
                   transcripts.hpp transcripts.cpp \
                   trie.hpp

isolator_LDADD = gtf/libgtf.la \
                 hat-trie/libhat-trie.la \
                 samtools/libbam.la \
                 sqlite3/libsqlite3.la \
                 -lz \
                 -lgsl \
                 -lgslcblas \
                 $(BOOST_THREAD_LIBS) \
                 $(BOOST_SYSTEM_LIBS) \
                 $(BOOST_TIMER_LIBS)

