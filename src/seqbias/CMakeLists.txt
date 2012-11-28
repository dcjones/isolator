
include_directories(${isolator_SOURCE_DIR}/src)

add_library(seqbias
            kmer_matrix.cpp
            motif.cpp
            sequencing_bias.cpp
            twobitseq.cpp
            pos_table.c)

