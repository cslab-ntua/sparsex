/*
 * prfcnt_xeon.h -- prfcnt Netburst implementation
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * Originally derived from gtsouk's code
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __PRFCNT_XEON_H__
#define __PRFCNT_XEON_H__

#include "prfcnt_common.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#include <linux/types.h> /* this is for __u64 and __u32 */

/*
 * ESCR registers fields
 */
#define ESCR_EVENT_SELECT_SHIFT  25
#define ESCR_EVENT_MASK_SHIFT    9
#define ESCR_TAG_VAL_SHIFT       5
#define ESCR_TAG_ENABLE_SHIFT    4
#define ESCR_OS_SHIFT            3
#define ESCR_USR_SHIFT           2

#define ESCR_T0_OS_SHIFT         3
#define ESCR_T0_USR_SHIFT        2
#define ESCR_T1_OS_SHIFT         1
#define ESCR_T1_USR_SHIFT        0

#define ESCR_SET(what, val) (val << (ESCR_ ## what ## _SHIFT) )

/*
 * CCCR registers fields
 */
#define CCCR_OVF_SHIFT           31
#define CCCR_CACSCADE_SHIFT      30
#define CCCR_OVF_PMI_SHIFT       26
#define CCCR_FORCE_OVF_SHIFT     25
#define CCCR_EDGE_SHIFT          24
#define CCCR_THRESHOLD_SHIFT     20
#define CCCR_COMPLEMENT_SHIFT    19
#define CCCR_COMPARE_SHIFT       18
#define CCCR_ESCR_SELECT_SHIFT   13
#define CCCR_ENABLE              12


#define CCCR_OVF_PMI_T0_SHIFT    26
#define CCCR_OVF_PMI_T1_SHIFT    27
#define CCCR_ACTIVE_THREAD_SHIFT 16

#define CCCR_SET(what, val) (val << (CCCR_ ## what ## _SHIFT) )
#define CCCR_EDGE_TRIGGER (CCCR_SET(COMPARE,1) | CCCR_SET(EDGE,1) | CCCR_SET(THRESHOLD,0))

/*
 * functions for handling counters:
 * --------------------------------
 *  prfcnt_cntr_setup: configure a counter for non-retirement events
 *  prfcnt_cntr_start: start counting on a counter
 *  prfcnt_cntr_pause: pause counting on a counter
 */

#define CCCR_ENABLE_BIT 12   /* ia32-3: 15.10.3, 15.11.2 (HT) */


#define IA32_PEBS_ENABLE      0x03f1
#define MSR_PEBS_MATRIX_VERT  0x03f2

struct prfcnt_xeon_event {
	__u32        cntr_nr;
	__u64        cntr_addr;
	__u64        escr_data, cccr_data;
	__u32        escr_addr, cccr_addr;

	char         *desc;
};
typedef struct prfcnt_xeon_event prfcnt_event_t;

struct prfcnt_xeon_handle {
	int             msr_fd;
	prfcnt_event_t  *events;
	unsigned int    cpu;
	unsigned long   flags;
};
typedef struct prfcnt_xeon_handle prfcnt_t;

enum {
	PRFCNT_AT_REP_L1_LMISS = 0x00, /* 1st Level cache load miss retired */
	PRFCNT_AT_REP_L2_LMISS,        /* 2nd level cache load miss retired */
	PRFCNT_AT_REP_DTLB_LMISS,      /* DTLB load miss retired */
	PRFCNT_AT_REP_DTLB_SMISS,      /* DTLB store miss retired */
	PRFCNT_AT_REP_DTLB_MISS        /* DTLB all miss retired */
	/* more ... */
};

enum {
	EVENT_L1_MISSES = 0,
	EVENT_L2_MISSES,
	EVENT_L2_ACCESSES,
	EVENT_IRET,
	EVENT_URET,
	EVENT_IOQ_ALLOC,
	//EVENT_IOQ_ENTRIES,
	EVENT_MCLEAR,
	EVENT_TC_MISSES,
	EVENT_MP_BRANCES,
	EVENT_END
};

static prfcnt_event_t __evnts[] = {
	{
		/*
		 * L1 MISSES
		 * escr: MSR_CRU_ESCR2
		 */
		.cntr_nr      = 12,

		.cntr_addr    = 0x030c,
		.cccr_addr    = 0x036c,
		.escr_addr    = 0x03cc,

		.escr_data    = ( ESCR_SET(EVENT_SELECT, 0x09) |
		                  ESCR_SET(EVENT_MASK,  (1<<0)) |
		                  ESCR_SET(TAG_ENABLE, 1) ),
		.cccr_data    = ( CCCR_SET(ESCR_SELECT, 0x05) |
		                  CCCR_SET(ACTIVE_THREAD, 3) ),

		.desc         = "l1 misses"
	},{
		/*
		 * L2 MISSES
		 * event:    BSQ_cache_reference
		 * escr:     BSU_CR_ECR0
		 */
		.cntr_nr    = 0,

		.cntr_addr  = 0x300,
		.escr_addr  = 0x3a0,
		.cccr_addr  = 0x360,

		/*
		 * escr event mask bits
		 */
		#define L2_RD_HITS (1<<0) /* shared    */
		#define L2_RD_HITE (1<<1) /* exclusive */
		#define L2_RD_HITM (1<<2) /* modified  */
		#define L3_RD_HITS (1<<3)
		#define L3_RD_HITE (1<<4)
		#define L3_RD_HITM (1<<5)
		#define L2_RD_MISS (1<<8)
		#define L3_RD_MISS (1<<9)
		#define L2_WR_MISS (1<<10)

		#define BSQ_ESCR_EVENT_MASK (L2_RD_MISS|L2_WR_MISS)
		#define BSQ_ESCR_EVENT_DESC "l2 misses"

		.escr_data  = ( ESCR_SET(EVENT_SELECT, 0x0c) |
		                ESCR_SET(EVENT_MASK, BSQ_ESCR_EVENT_MASK) ),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT,0x07)  |
		                CCCR_SET(ACTIVE_THREAD,3)),

		.desc       = BSQ_ESCR_EVENT_DESC

		#undef BSQ_ESCR_EVENT_MASK
		#undef BSQ_ESCR_EVENT_DESC
	},{
		/*
		 * L2 ACCESSES
		 * event:    BSQ_cache_reference
		 * escr:     BSU_CR_ECR1
		 */
		.cntr_nr    = 2,

		.cntr_addr  = 0x302,
		.escr_addr  = 0x3a1,
		.cccr_addr  = 0x362,

		/*
		 * escr event mask bits
		 */
		#define L2_RD_HITS (1<<0) /* shared    */
		#define L2_RD_HITE (1<<1) /* exclusive */
		#define L2_RD_HITM (1<<2) /* modified  */
		#define L3_RD_HITS (1<<3)
		#define L3_RD_HITE (1<<4)
		#define L3_RD_HITM (1<<5)
		#define L2_RD_MISS (1<<8)
		#define L3_RD_MISS (1<<9)
		#define L2_WR_MISS (1<<10)

		#define BSQ_ESCR_EVENT_MASK (L2_RD_MISS|L2_WR_MISS|L2_RD_HITS|L2_RD_HITE|L2_RD_HITM)
		#define BSQ_ESCR_EVENT_DESC "l2 accesses"

		.escr_data  = ( ESCR_SET(EVENT_SELECT, 0x0c) |
		                ESCR_SET(EVENT_MASK, BSQ_ESCR_EVENT_MASK) ),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT,0x07)  |
		                CCCR_SET(ACTIVE_THREAD,3)),

		.desc       = BSQ_ESCR_EVENT_DESC
	},{
		/*
		 * INSTRUCTIONS RETIRED
		 * escr: MSR_CRU_ESCR0
		 */
		.cntr_nr    = 13,

		.cntr_addr  = 0x030d,
		.escr_addr  = 0x03b8,
		.cccr_addr  = 0x036d,

		#define INSTR_NBOGUS_NTAG (1<<0) /* non-bogus && non-tagged */
		#define INSTR_NBOGUS_TAG  (1<<1) /* non-bogus && tagged */
		#define INSTR_BOGUS_NTAG  (1<<2) /* bogus && non-tagged */
		#define INSTR_BOGUS_TAG   (1<<3) /* bogus && tagged */
		#define INSTR_NBOGUS      (INSTR_NBOGUS_NTAG|INSTR_NBOGUS_TAG)
		#define INSTR_BOGUS       (INSTR_BOGUS_NTAG|INSTR_BOGUS_TAG)

		#define INSTR_ESCR_EVENT_MASK INSTR_NBOGUS
		#define INSTR_ESCR_EVENT_DESC "iret (nbogus)"


		.escr_data  = ( ESCR_SET(EVENT_SELECT,0x02) |
		                ESCR_SET(EVENT_MASK, INSTR_ESCR_EVENT_MASK ) ),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT,0x4) |
		                CCCR_SET(ACTIVE_THREAD, 3) ),

		.desc       = INSTR_ESCR_EVENT_DESC

	},{
		/*
		 * UOPS RETIRED
		 * escr: MSR_CRU_ESCR1
		 */
		.cntr_nr    = 14,

		.cntr_addr  = 0x030e,
		.escr_addr  = 0x03b9,
		.cccr_addr  = 0x036e,

		#define UOPS_NBOGUS (1<<0)
		#define UOPS_BOGUS  (1<<1)

		#define UOPS_ESCR_EVENT_MASK UOPS_NBOGUS
		#define UOPS_ESCR_EVENT_DESC "uops (nbogus)"

		.escr_data  = ( ESCR_SET(EVENT_SELECT,0x01) |
		                ESCR_SET(EVENT_MASK, UOPS_ESCR_EVENT_MASK)),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT,0x4) |
		                CCCR_SET(ACTIVE_THREAD, 3) ),

		.desc       = UOPS_ESCR_EVENT_DESC
	},{
		/*
		 * IOQ ALLOCATION
		 * escr: MSR_FSB_ESCR0
		 */
		.cntr_nr    = 1,

		.cntr_addr  = 0x0301,
		.escr_addr  = 0x03a2,
		.cccr_addr  = 0x0361,

		#define IOQ_ALL_RD   (1<<5)
		#define IOQ_ALL_WR   (1<<6)
		#define IOQ_MEM_UC   (1<<7)
		#define IOQ_MEM_WC   (1<<8)
		#define IOQ_MEM_WT   (1<<9)
		#define IOQ_MEM_WP   (1<<10)
		#define IOQ_MEM_WB   (1<<11)
		#define IOQ_MEM_ALL  (IOQ_MEM_UC|IOQ_MEM_WC|IOQ_MEM_WT|IOQ_MEM_WP|IOQ_MEM_WB)
		#define IOQ_OWN      (1<<13)
		#define IOQ_OTHER    (1<<14)
		#define IOQ_PREFETCH (1<<15)

		#define IOQ_ESCR_EVENT_MASK (IOQ_ALL_RD|IOQ_MEM_ALL|IOQ_OWN)
		#define IOQ_ESCR_EVENT_DESC "(R)"



		.escr_data = ( ESCR_SET(EVENT_SELECT,0x03) |
		               ESCR_SET(EVENT_MASK, IOQ_ESCR_EVENT_MASK)),
		.cccr_data = ( CCCR_SET(ESCR_SELECT, 0x6) |
		               CCCR_SET(ACTIVE_THREAD, 3) |
		               CCCR_EDGE_TRIGGER ),

		.desc       = "IOQ_alloc " IOQ_ESCR_EVENT_DESC
	},
	#if 0
	{
		/*
		 * IOQ ENTRIES
		 * escr: MSR_FSB_ESCR1
		 */
		.cntr_nr    = 2,

		.cntr_addr  = 0x0302,
		.escr_addr  = 0x03a3,
		.cccr_addr  = 0x0362,

		.escr_data  = ( ESCR_SET(EVENT_SELECT,0x1a) |
		                ESCR_SET(EVENT_MASK, IOQ_ESCR_EVENT_MASK )),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT, 0x6) |
		                CCCR_SET(ACTIVE_THREAD, 3) ),

		.desc       = "IOQ_entries " IOQ_ESCR_EVENT_DESC
	},
	#endif
	{
		/*
		 * MACHINE CLEAR
		 * escr: MSR_CRU_ESCR3
		 */
		.cntr_nr    = 15,

		#define MCLEAR_CLEAR    (1<<0)
		#define MCLEAR_MOCLEAR  (1<<2) /* memory ordering */
		#define MCLEAR_SMCLEAR  (1<<6) /* self modifying */

		#define MCLEAR_ESCR_EVENT_MASK (MCLEAR_MOCLEAR)
		#define MCLEAR_ESCR_EVENT_DESC  "(mo)"

		.cntr_addr  = 0x030f,
		.escr_addr  = 0x03cd,
		.cccr_addr  = 0x036f,

		.escr_data  = ( ESCR_SET(EVENT_SELECT,0x02) |
		                ESCR_SET(EVENT_MASK, MCLEAR_ESCR_EVENT_MASK)),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT, 0x05) |
		                CCCR_SET(ACTIVE_THREAD, 3) ),

		.desc       = "mclear " MCLEAR_ESCR_EVENT_DESC
	},{
		/*
		 * TC Misses
		 * escr: MSR_BPU_ESCR1
		 */
		.cntr_nr    = 3,

		.cntr_addr  = 0x0303,
		.escr_addr  = 0x03b3,
		.cccr_addr  = 0x0363,


		.escr_data  = ( ESCR_SET(EVENT_SELECT, 0x18) |
		                ESCR_SET(EVENT_MASK, 1)),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT, 0x00) |
		                CCCR_SET(ACTIVE_THREAD, 3) ),

		.desc       = "TC misses"
	},{
		/*
		 * MP Branches
		 * escr: MSR_TBPU_ESCR0
		 */
		.cntr_nr    = 4,

		.cntr_addr  = 0x0304,
		.escr_addr  = 0x03c2,
		.cccr_addr  = 0x0364,

		#define MP_BRANCH_COND  (1<<1) /* conditional */
		#define MP_BRANCH_CALL  (1<<2) /* indirect call */
		#define MP_BRANCH_RET   (1<<3) /* return branch */
		#define MP_BRANCH_IND   (1<<4) /* return, indirect call/jump */

		#define MP_BRANCH_ESCR_EVENT_MASK (MP_BRANCH_COND|MP_BRANCH_IND)
		#define MP_BRANCH_ESCR_EVENT_DESC "(ind+cond)"

		.escr_data  = ( ESCR_SET(EVENT_SELECT, 0x05) |
		                ESCR_SET(EVENT_MASK, MP_BRANCH_ESCR_EVENT_MASK)),
		.cccr_data  = ( CCCR_SET(ESCR_SELECT, 0x02) |
		                CCCR_SET(ACTIVE_THREAD, 3) ),

		.desc       = "MP Branches " MP_BRANCH_ESCR_EVENT_DESC
	}
};

/*
 * fields of ia32_pebs_enable msr which are related to at
 * retirements replay event.
 *
 *   [0..12]: replay event metric,
 *   [24]   : enable replay tagging
 *  HT:
 *   [25]   : enable pebs for logical processor my thread runs
 *   [26]   : enable pebs for the other logical processor
 *  NO HT:
 *   [25]   : enables pebs
 *   [26]   : reserverd
 */
static const struct {
	__u64 ia32_pebs_enable_data;
	__u64 msr_pebs_matrix_vert_data;
} prfcnt_rep_metrics[] = {
	{	/* PRFCNT_AT_REP_L1_LMISS */
		.ia32_pebs_enable_data     = (1<<0) | (1<<24),  (1<<25),
		.msr_pebs_matrix_vert_data = (1<<0)
	},{	/* PRFCNT_AT_REP_L2_LMISS */
		.ia32_pebs_enable_data     = (1<<1) | (1<<24) | (1<<25),
		.msr_pebs_matrix_vert_data = (1<<0)
	}
	/* more ... */
};

static inline void prfcnt_cntr_setup(prfcnt_t *handle, prfcnt_event_t *event)
{
	event->cccr_data &= ~(1<<CCCR_ENABLE_BIT);

	prfcnt_wrmsr(handle->msr_fd, event->cccr_addr, event->cccr_data);
	prfcnt_wrmsr(handle->msr_fd, event->escr_addr, event->escr_data);
	prfcnt_wrmsr(handle->msr_fd, event->cntr_addr, 0);

}

static inline void prfcnt_cntr_at_rep_setup(prfcnt_t *handle,
                                            prfcnt_event_t *event,
                                            unsigned int rep_metric)
{
	prfcnt_wrmsr(handle->msr_fd, IA32_PEBS_ENABLE,
	             prfcnt_rep_metrics[rep_metric].ia32_pebs_enable_data);
	prfcnt_wrmsr(handle->msr_fd, MSR_PEBS_MATRIX_VERT,
	             prfcnt_rep_metrics[rep_metric].msr_pebs_matrix_vert_data);

	prfcnt_cntr_setup(handle, event);
}

static inline void prfcnt_cntr_start(prfcnt_t *handle, prfcnt_event_t *event)
{
	event->cccr_data |= (1<<CCCR_ENABLE_BIT);
	prfcnt_wrmsr(handle->msr_fd,event->cccr_addr,event->cccr_data);
}

static inline void prfcnt_cntr_pause(prfcnt_t *handle, prfcnt_event_t *event)
{
	event->cccr_data &= ~(1<<CCCR_ENABLE_BIT);
	prfcnt_wrmsr(handle->msr_fd,event->cccr_addr,event->cccr_data);
}


#define prfcnt_cntr_readpmc(i) \
	prfcnt_rdpmc(__evnts[i].cntr_nr)

#define PRFCNT_FL_T0 (1<<0)  /* count events on T0 */
#define PRFCNT_FL_T1 (1<<1)  /* count events on T1 */

static inline void prfcnt_cntr_set_t0(prfcnt_event_t *event)
{
	event->escr_data |= ESCR_SET(T0_USR,1);
}

static inline void prfcnt_cntr_set_t1(prfcnt_event_t *event)
{
	event->escr_data |= ESCR_SET(T1_USR,1);
}

static inline void prfcnt_init(prfcnt_t *handle, int cpu, unsigned long flags)
{
	int i;

	handle->cpu    = cpu;
	handle->flags  = flags;
	handle->events = __evnts;

	__prfcnt_init(cpu, &handle->msr_fd);

	for ( i=0; i<EVENT_END; i++){
		prfcnt_event_t *event = &handle->events[i];

		if (flags & PRFCNT_FL_T0)
			prfcnt_cntr_set_t0(event);

		if (flags & PRFCNT_FL_T1)
			prfcnt_cntr_set_t1(event);

		switch (i) {
			case EVENT_L1_MISSES:
			prfcnt_cntr_at_rep_setup(handle, event,
			                         PRFCNT_AT_REP_L1_LMISS);
			break;

			default:
			prfcnt_cntr_setup(handle,event);
		}
	}
}

static inline void prfcnt_shut(prfcnt_t *handle)
{
	__prfcnt_shut(handle->msr_fd);
}

static inline void prfcnt_start(prfcnt_t *handle)
{
	int i;

	for (i=0; i<EVENT_END; i++){
		prfcnt_cntr_start(handle,&handle->events[i]);
	}
}

static inline void prfcnt_pause(prfcnt_t *handle)
{
	int i;

	for (i=0; i<EVENT_END; i++){
		prfcnt_cntr_pause(handle,&handle->events[i]);
	}
}

static inline __u64  *prfcnt_readstats_rdpmc(prfcnt_t *handle, __u64 *stats)
{
	int i;

	for (i=0; i < EVENT_END; i++){
		stats[i] = prfcnt_rdpmc(handle->events[i].cntr_nr);
	}

	return stats;
}

static inline __u64  *prfcnt_readstats_rdmsr(prfcnt_t *handle, __u64 *stats)
{
	int i;

	for (i=0; i < EVENT_END; i++){
		prfcnt_rdmsr(handle->msr_fd,handle->events[i].cntr_addr, &stats[i]);
	}

	return stats;
}

#define prfcnt_readstats prfcnt_readstats_rdmsr

static inline void prfcnt_report(prfcnt_t *handle)
{
	__u64 stats[EVENT_END];
	int   i;

	prfcnt_readstats(handle,stats);

	for (i=0; i < EVENT_END; i++){
		printf("\n\t%25s:%16llu", handle->events[i].desc, stats[i]);
	}
	printf("\n");
}

#endif /* __PRFCNT_H__ */
