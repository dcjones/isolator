
SUBDIRS = gtf hat-trie samtools

bin_PROGRAMS = isolator

isolator_SOURCES = isolator.cpp \
                   common.hpp common.cpp \
                   transcripts.hpp transcripts.cpp \
                   trie.hpp

isolator_LDADD = gtf/libgtf.la \
                 hat-trie/libhat-trie.la \
                 samtools/libbam.la \
                 -lz \
                 -lgsl \
                 -lgslcblas \
                 $(BOOST_THREAD_LIB)
