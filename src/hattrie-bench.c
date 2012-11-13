
#include <stdio.h>

#include "samtools/sam.h"
#include "hat-trie/hat-trie.h"

int main(int argc, char* argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: hattrie-bench reads.bam\n");
        return 1;
    }

    samfile_t* f = samopen(argv[1], "rb", NULL);
    if (f == NULL) {
        fprintf(stderr, "Can't open BAM file %s.", argv[1]);
        return 1;
    }
    bam1_t* b = bam_init1();
    hattrie_t* t = hattrie_create();
    while (samread(f, b) > 0) {
        value_t* v = hattrie_get(t, bam1_qname(b), strlen(bam1_qname(b)));
        *v = 1;
    }

    hattrie_free(t);
    bam_destroy1(b);
    samclose(f);
    return 0;
}


