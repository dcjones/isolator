
SUBDIRS = hat-trie samtools

bin_PROGRAMS = isolator

isolator_SOURCES = isolator.cpp \
                   transcripts.hpp transcripts.cpp \
                   trie.hpp

isolator_LDADD = hat-trie/libhat-trie.la \
				 samtools/libbam.la \
                 -lz \
                 -lgsl \
                 -lgslcblas \
                 $(BOOST_THREAD_LIB)
