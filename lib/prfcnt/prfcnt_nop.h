#ifndef __PRFCNT_NOP_H__
#define __PRFCNT_NOP_H__


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif


#define PRFCNT_FL_T0 0
#define PRFCNT_FL_T1 0

#include <stdlib.h>
#include <inttypes.h>

typedef int prfcnt_t;
typedef int prfcnt_event_t;

static inline void prfcnt_cntr_setup(prfcnt_t *handle, prfcnt_event_t *event)
{
}

static inline void prfcnt_cntr_start(prfcnt_t *handle, prfcnt_event_t *event)
{
}

static inline void prfcnt_cntr_pause(prfcnt_t *handle, prfcnt_event_t *event)
{
}


static inline void prfcnt_init(prfcnt_t *handle, int cpu, unsigned long flags)
{
}

static inline void prfcnt_shut(prfcnt_t *handle)
{
}

static inline void prfcnt_start(prfcnt_t *handle)
{
}

static inline void prfcnt_pause(prfcnt_t *handle)
{
}

static inline uint64_t *prfcnt_readstats(prfcnt_t *handle, uint64_t *stats)
{
	return NULL;
}

static inline void prfcnt_report(prfcnt_t *handle)
{
}

#endif /* __PRFCNT_NOP_H__ */
