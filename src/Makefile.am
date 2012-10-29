
SUBDIRS = gtf hat-trie samtools sqlite3

bin_PROGRAMS = isolator

isolator_SOURCES = isolator.cpp \
                   common.hpp common.cpp \
                   constants.hpp constants.cpp \
                   logger.hpp logger.cpp \
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

