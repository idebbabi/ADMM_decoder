#ifndef rdtsc_timer_imen
#define rdtsc_timer_imen

static __inline__ unsigned long long timer() {
    unsigned hi, lo;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long) lo) | (((unsigned long long) hi) << 32);
}

#endif