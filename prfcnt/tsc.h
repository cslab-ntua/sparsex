#ifndef __TSC__
#define __TSC__

#include <stdio.h>

#include <linux/types.h>

struct tsc {
	__u64	ticks;
	__u64	last;
}; 
typedef struct tsc tsc_t;

#if defined(__i386__) || defined(__x86_64__)
static inline __u64 get_ticks()
{
	__u32 hi,low;
	__u64 ret;

	__asm__ __volatile__ ("rdtsc" : "=a"(low), "=d"(hi));

	ret = hi;
	ret <<= 32;
	ret |= low;
	
	return ret;
}
#elif defined(__ia64__)
#include <asm/intrinsics.h>
static inline __u64 get_ticks()
{
	__u64 ret = ia64_getreg(_IA64_REG_AR_ITC);	
	ia64_barrier();

	return ret;
}
#else 
#error "dont know how to count ticks"
#endif

static inline void tsc_init(tsc_t *tsc)
{
	tsc->ticks = 0;
	tsc->last  = 0;
}

static inline void tsc_shut(tsc_t *tsc)
{
}

static inline void tsc_start(tsc_t *tsc)
{
	tsc->last = get_ticks();
}

static inline void tsc_pause(tsc_t *tsc)
{
	tsc->ticks += (get_ticks() - tsc->last);
}


#ifdef CPU_MHZ
static inline double __tsc_getsecs(__u64 ticks)
{
	return (ticks/(1000000*CPU_MHZ));
}
static inline double tsc_getsecs(tsc_t *tsc)
{
	return __tsc_getsecs(tsc->ticks);
}
#else
#warning no CPU_MHZ defined 
#endif

static inline __u64 tsc_getticks(tsc_t *tsc)
{
	return tsc->ticks;
}

static inline void tsc_report(tsc_t *tsc)
{
	__u64 ticks = tsc->ticks;
	
	printf("ticks : %llu\n", (long long unsigned)ticks);
	#ifdef CPU_MHZ
	printf("time  : %lf (sec)\n", __tsc_getsecs(ticks));
	#endif
}

#endif
