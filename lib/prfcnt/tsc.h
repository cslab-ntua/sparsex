#ifndef __TSC__
#define __TSC__

/*
 * Code to count ticks
 * Kornilios Kourtis <kkourt@cslab.ece.ntua.gr>
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <inttypes.h>

struct tsc {
	uint64_t	ticks;
	uint64_t	last;
};
typedef struct tsc tsc_t;

#if defined(__i386__) || defined(__x86_64__)
static inline uint64_t get_ticks()
{
	uint32_t hi,low;
	uint64_t ret;

	__asm__ __volatile__ ("rdtsc" : "=a"(low), "=d"(hi));

	ret = hi;
	ret <<= 32;
	ret |= low;

	return ret;
}
#elif defined(__ia64__)
#include <asm/intrinsics.h>
static inline uint64_t get_ticks()
{
	uint64_t ret = ia64_getreg(_IA64_REG_AR_ITC);
	ia64_barrier();

	return ret;
}
#elif defined(__sparc__)
// linux-2.6.28/arch/sparc64/kernel/time.c
static inline uint64_t get_ticks()
{
	uint64_t t;
	__asm__ __volatile__ (
		"rd     %%tick, %0\n\t"
		"mov    %0,     %0"
		: "=r"(t)
	);
	return t & (~(1UL << 63));
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

/*
 * This is ugly, but we need a way to determine the CPU's frequency, so that we
 * can translate ticks to seconds. It is better to do it at run-time to avoid
 * stupid errors (yes, it used to be set at compile time, and I wasted several
 * hours chasing the `bug').
 */
static double __getMhz(void)
{
#ifdef CPU_MHZ_SH
	double mhz;
	FILE *script;
	char buff[512], *endptr;
	int ret;

	script = popen(CPU_MHZ_SH, "r");
	if (script == NULL){
			perror("popen");
			goto error;
	}

	ret = fread(buff, 1, sizeof(buff), script);
	if (!ret){
		perror("fread");
		goto error;
	}

	mhz = strtod(buff, &endptr);
	if (endptr == buff){
		fprintf(stderr, "strtod failed\n");
		goto error;
	}

	return mhz;

error:
#endif
	fprintf(stderr,
	        "%s@%s:%d: ERROR: Dont now how to calculate MhZ.\n"
	        "Please set CPU_MHZ_SH to a script that returns MHz when compiling\n"
	        "Aborting...\n",
	         __FUNCTION__, __FILE__,__LINE__);
	exit(1);
	return 0.0;
}

static inline double getMhz(void)
{
	static double mhz = 0.0;
	if (mhz == 0.0){
		mhz = __getMhz();
	}
	return mhz;
}

static inline double __tsc_getsecs(uint64_t ticks)
{
	return (ticks/(1000000.0*getMhz()));
}
static inline double tsc_getsecs(tsc_t *tsc)
{
	return __tsc_getsecs(tsc->ticks);
}

static uint64_t tsc_getticks(tsc_t *tsc)
{
	return tsc->ticks;
}

static inline void tsc_spinticks(uint64_t ticks)
{
	uint64_t t0;
	//uint64_t spins = 0;
	t0 = get_ticks();
	for (;;){
		if (get_ticks() - t0 >= ticks)
			break;
		//spins++;
	}
	//printf("spins=%lu\n", spins);
}

static inline void tsc_report(tsc_t *tsc)
{
	uint64_t ticks = tsc_getticks(tsc);

	printf("ticks : %llu\n", (long long unsigned)ticks);
	printf("time  : %lf (sec)\n", __tsc_getsecs(ticks));
}

#endif
