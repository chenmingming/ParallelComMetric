#ifndef RDTSC_H_
#define RDTSC_H_

inline unsigned long long rdtsc(void)
{
        unsigned long long result=0;
        unsigned a, d;

        do {
                __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d));
                result = ((unsigned long long)a) | (((unsigned long long)d) << 32);
        } while (__builtin_expect ((int) result == -1, 0));
        return result;
}

#endif
