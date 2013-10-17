
#include <cstdio>


static void cpuid(unsigned int op, unsigned int* eax,
                  unsigned int* ecx, unsigned int* edx)
{
     __asm__("pushq %%rbx\n\tcpuid\n\tpopq %%rbx"
         : "=a" (*eax), "=c" (*ecx), "=d" (*edx)
         : "a" (op));
}


static void xgetbv(unsigned int op, unsigned int* eax, unsigned int* edx)
{
    __asm__(".byte 0x0f, 0x01, 0xd0"
        : "=a" (*eax), "=d" (*edx)
        : "c" (op));
}


bool cpu_has_sse2()
{
    unsigned int eax, ecx, edx;
    cpuid(0, &eax, &ecx, &edx);
    if (eax == 0) {
        return false;
    }

    cpuid(1, &eax, &ecx, &edx);
    return edx & 0x04000000;
}


bool cpu_has_sse4()
{
    unsigned int eax, ecx, edx;
    cpuid(0, &eax, &ecx, &edx);
    if (eax == 0) {
        return false;
    }

    cpuid(1, &eax, &ecx, &edx);
    return ecx & 0x00080000;
}


bool cpu_has_avx()
{
    unsigned int eax, ecx, edx;
    cpuid(0, &eax, &ecx, &edx);
    if (eax == 0) {
        return false;
    }

    cpuid(1, &eax, &ecx, &edx);
    if ((ecx & 0x18000000) == 0x18000000 ) {
        xgetbv(0, &eax, &edx);
        return eax & 0x6;
    }
    else return false;
}


#if 0
int main()
{
    unsigned int eax, ecx, edx;
    cpuid(0, &eax, &ecx, &edx);

    printf("cpu features: ");
    if (cpu_has_sse2()) printf(" sse2");
    if (cpu_has_sse4()) printf(" sse4");
    if (cpu_has_avx()) printf(" avx");
    printf("\n");

    return 0;
}
#endif


