
#include <cctype>
#include <cstdio>
#include <cstring>

#include "common.hpp"


strand_t other_strand(strand_t s)
{
    if      (s == strand_na)  return strand_na;
    else if (s == strand_pos) return strand_neg;
    else                      return strand_pos;
}


static int seqname_num(const char* u)
{
    int num = 0;
    while(*u != '\0') {
        if (isdigit(*u)) {
            sscanf(u, "%d", &num);
            break;
        }

        ++u;
    }

    return num;
}


int seqname_compare(const char* u, const char* v)
{
    int x = seqname_num(u);
    int y = seqname_num(v);

    if (x != y) return x - y;
    else        return strcmp(u, v);
}

