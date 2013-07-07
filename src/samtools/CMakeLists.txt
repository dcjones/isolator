
add_library(samtools
            bam.c
            bam_aux.c
            bam_import.c
            bam_index.c
            bam_pileup.c
            bgzf.c
            faidx.c
            knetfile.c
            kstring.c
            razf.c
            sam.c
            sam_header.c
            samtools_extra.c)

add_definitions(-D_USE_KNETFILE)

# Not my code. Supress compiler warnings.
set_target_properties(samtools PROPERTIES COMPILE_FLAGS "-w")

