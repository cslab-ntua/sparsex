#ifndef __PRFCNT_SIMPLE_H__
#define __PRFCNT_SIMPLE_H__

/*
 * prfcnt_simple.h
 *
 * Common functions for Core and AMD Architecture
 */

static inline void prfcnt_cntr_setup(prfcnt_t *handle, prfcnt_event_t *event)
{
	__u64        evntsel_data = event->evntsel_data;
	unsigned int counter_addr = event->cntr_nr + PRFCNT_CNTR_BADDR; 
	unsigned int evntsel_addr = event->cntr_nr + PRFCNT_EVNT_BADDR;

	evntsel_data &= ~(1<<EVNTSEL_ENABLE_CNTR_SHIFT);

	prfcnt_wrmsr(handle->msr_fd, counter_addr, 0);
	prfcnt_wrmsr(handle->msr_fd, evntsel_addr, evntsel_data);
}

static inline void prfcnt_cntr_start(prfcnt_t *handle, prfcnt_event_t *event)
{
	__u64        evntsel_data = event->evntsel_data;
	unsigned int evntsel_addr = event->cntr_nr + PRFCNT_EVNT_BADDR;

	evntsel_data |= (1<<EVNTSEL_ENABLE_CNTR_SHIFT);
	prfcnt_wrmsr(handle->msr_fd, evntsel_addr, evntsel_data);
}

static inline void prfcnt_cntr_pause(prfcnt_t *handle, prfcnt_event_t *event)
{
	__u64        evntsel_data = event->evntsel_data;
	unsigned int evntsel_addr = event->cntr_nr + PRFCNT_EVNT_BADDR;

	evntsel_data &= ~(1<<EVNTSEL_ENABLE_CNTR_SHIFT);
	prfcnt_wrmsr(handle->msr_fd, evntsel_addr, evntsel_data);
}


static inline void prfcnt_init(prfcnt_t *handle, int cpu, unsigned long flags)
{
	int i;

	handle->cpu       = cpu;
	handle->flags     = flags;
	handle->events_nr = sizeof(__evnts_selected)/sizeof(int);
	handle->events    = malloc(handle->events_nr*sizeof(prfcnt_event_t *));
	if ( !handle->events ){
		perror("prfcnt_init");
		exit(1);
	}

	for (i=0; i < handle->events_nr ; i++){
		prfcnt_event_t *event;
		event             = &__evnts[__evnts_selected[i]];
		event->cntr_nr    = i;
		handle->events[i] = event;
	}

	__prfcnt_init(cpu, &handle->msr_fd);

	for (i=0; i < handle->events_nr; i++){
		prfcnt_cntr_setup(handle, handle->events[i]);
	}
}

static inline void prfcnt_shut(prfcnt_t *handle)
{
	__prfcnt_shut(handle->msr_fd);
}

static inline void prfcnt_start(prfcnt_t *handle)
{
	int i;

	for (i=0; i < handle->events_nr; i++){
		prfcnt_cntr_start(handle, handle->events[i]);
	}
}

static inline void prfcnt_pause(prfcnt_t *handle)
{
	int i;

	for (i=0; i < handle->events_nr; i++){
		prfcnt_cntr_pause(handle, handle->events[i]);
	}
}

static inline __u64 *prfcnt_readstats_rdmsr(prfcnt_t *handle, __u64 *stats)
{
	int i;
	unsigned int cntr_addr;

	for (i=0; i < handle->events_nr; i++){
		cntr_addr = handle->events[i]->cntr_nr + PRFCNT_CNTR_BADDR;
		prfcnt_rdmsr(handle->msr_fd, cntr_addr, stats + i);
	}

	return stats;
}

static inline __u64  *prfcnt_readstats_rdpmc(prfcnt_t *handle, __u64 *stats)
{
	int i;
	
	for (i=0; i < handle->events_nr; i++){
		stats[i] = prfcnt_rdpmc(handle->events[i]->cntr_nr);
	}

	return stats;
}

#define prfcnt_readstats prfcnt_readstats_rdmsr

static inline void prfcnt_report(prfcnt_t *handle)
{
	__u64 stats[EVENT_END];
	int   i;

	prfcnt_readstats(handle, stats);

	printf(" cpu %d:\n", handle->cpu);
	for (i=0; i < handle->events_nr; i++){
		printf("\t%50s:%12llu\n", handle->events[i]->desc, stats[i]);
	}
}

#endif /* __PRFCNT_SIMPLE_H__ */
