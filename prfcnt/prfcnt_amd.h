#ifndef __PRFCNT_AMD__
#define __PRFCNT_AMD__


#include "prfcnt_common.h"

#define PRFCNT_EVNT_BADDR 0xc0010000 /* event selectors base address */
#define PRFCNT_CNTR_BADDR 0xc0010004 /* perfomance counters base address */ 

/*
 * selector registers fields
 */
#define EVNTSEL_CNTMASK_SHIFT      24 
#define EVNTSEL_INVERT_SHIFT       23
#define EVNTSEL_ENABLE_CNTR_SHIFT  22
#define EVNTSEL_ENABLE_APIC_SHIFT  20
#define EVNTSEL_PIN_CTRL_SHIFT     19
#define EVNTSEL_EDGE_DETECT_SHIFT  18
#define EVNTSEL_OS_SHIFT           17
#define EVNTSEL_USR_SHIFT          16
#define EVNTSEL_UNIT_SHIFT         8
#define EVNTSEL_EVENT_SHIFT        0 

#define EVNTSEL_SET(what, val) ( val << (EVNTSEL_ ## what ## _SHIFT) )

#define EVNT_UNIT_MOESI_INVALID     (1<<0)
#define EVNT_UNIT_MOESI_SHARED      (1<<1)
#define EVNT_UNIT_MOESI_EXCLUSIVE   (1<<2)
#define EVNT_UNIT_MOESI_OWNER       (1<<3)
#define EVNT_UNIT_MOESI_MODIFIED    (1<<4)
#define EVNT_UNIT_MOESI_NOT_INVALID ( EVNT_UNIT_MOESI_SHARED    |\
                                      EVNT_UNIT_MOESI_EXCLUSIVE |\
                                      EVNT_UNIT_MOESI_OWNER     |\
                                      EVNT_UNIT_MOESI_MODIFIED )
#define EVNT_UNIT_MOESI_ALL         ( EVNT_UNIT_MOESI_INVALID | EVNT_UNIT_MOESI_NOT_INVALID) 
                                  
#define EVNT_UNIT_TLB_RELOAD (1<<2)
                                 
                                
#define EVNT_EVENT_DC_ACCESS     0x40
#define EVNT_EVENT_DC_MISS       0x41
#define EVNT_EVENT_DC_REFILL_L2  0x42
#define EVNT_EVENT_DC_REFILL_SYS 0x43
#define EVNT_EVENT_DC_COPYBACK   0x44
#define EVNT_EVENT_DTLB1_MISS    0x45
#define EVNT_EVENT_DTLB2_MISS    0x46
#define EVNT_EVENT_INT_L2_REQ    0x7d
#define EVNT_EVENT_INT_L2_MISS   0x7e
#define EVNT_EVENT_UOPS          0xc1
#define EVNT_EVENT_IC_REFILL_L2  0x82
#define EVNT_EVENT_IC_REFILL_SYS 0x83


struct prfcnt_amd_event {
	__u32   cntr_nr;
	__u64   evntsel_data;
	char    *desc;
};
typedef struct prfcnt_amd_event prfcnt_event_t;

struct prfcnt_amd_handle { 
	int                msr_fd;
	prfcnt_event_t     **events;
	unsigned int       events_nr;
	unsigned int       cpu;
	unsigned long      flags;
};
typedef struct prfcnt_amd_handle prfcnt_t;

enum {
	EVENT_DC_REFILLS_L2_ALL = 0,
	EVENT_DC_REFILLS_L2_INVALID,
	EVENT_DC_REFILLS_L2_SHARED,
	EVENT_DC_REFILLS_L2_EXCLUSIVE,
	EVENT_DC_REFILLS_L2_NOT_INVALID,
	EVENT_DC_REFILLS_SYS_ALL,
	EVENT_IC_REFILLS_L2_ALL,
	EVENT_IC_REFILLS_SYS_ALL,
	EVENT_DC_REFILLS_SYS_INVALID,
	EVENT_DC_REFILLS_SYS_SHARED,
	EVENT_DC_REFILLS_SYS_EXCLUSIVE,
	EVENT_DC_REFILLS_SYS_OWNER,
	EVENT_DC_REFILLS_SYS_MODIFIED,
	EVENT_L2_REQS_TLB_RELOAD,
	EVENT_L2_MISSES_TLB_RELOAD,
	EVENT_END
};

static prfcnt_event_t __evnts[] = {
	[EVENT_DC_REFILLS_L2_ALL] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_L2) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_ALL)      |
		                 EVNTSEL_SET(USR, 1)                         |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills L2 (all)"
	},
	[EVENT_DC_REFILLS_L2_INVALID] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_L2) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_INVALID)  |
		                 EVNTSEL_SET(USR, 1)                         |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills L2 (invalid)"
	},
	[EVENT_DC_REFILLS_L2_SHARED] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_L2) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_SHARED)  |
		                 EVNTSEL_SET(USR, 1)                         |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills L2 (shared)"
	},
	[EVENT_DC_REFILLS_L2_EXCLUSIVE] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_L2) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_EXCLUSIVE)  |
		                 EVNTSEL_SET(USR, 1)                         |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills L2 (exclusive)"
	},
	[EVENT_DC_REFILLS_L2_NOT_INVALID] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_L2) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_NOT_INVALID)  |
		                 EVNTSEL_SET(USR, 1)                         |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills L2 (not invalid)"
	},
	
	[EVENT_DC_REFILLS_SYS_ALL] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_SYS) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_ALL)     |
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills Sys (all)"
	},

	[EVENT_IC_REFILLS_L2_ALL] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_IC_REFILL_L2) | 
		                 EVNTSEL_SET(USR, 1)                         |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "IC refills L2 (all)"
	},
	
	[EVENT_IC_REFILLS_SYS_ALL] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_IC_REFILL_SYS) | 
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "IC refills Sys (all)"
	},
	
	[EVENT_DC_REFILLS_SYS_INVALID] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_SYS) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_INVALID) |
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills Sys (invalid)"
	},

	[EVENT_DC_REFILLS_SYS_SHARED] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_SYS) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_SHARED)  |
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills Sys (shared)"
	},

	[EVENT_DC_REFILLS_SYS_EXCLUSIVE] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_SYS)   | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_EXCLUSIVE) |
		                 EVNTSEL_SET(USR, 1)                          |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills Sys (exclusive)"
	},
	
	[EVENT_DC_REFILLS_SYS_OWNER] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_SYS) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_OWNER)   |
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills Sys (owner)"
	},

	[EVENT_DC_REFILLS_SYS_MODIFIED] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_DC_REFILL_SYS) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_MOESI_MODIFIED) |
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "DC refills Sys (modified)"
	},
	[EVENT_L2_REQS_TLB_RELOAD] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_INT_L2_REQ) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_TLB_RELOAD) |
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "L2 internal reqs (TLB)"
	},
	[EVENT_L2_MISSES_TLB_RELOAD] = {
		.evntsel_data  = EVNTSEL_SET(EVENT, EVNT_EVENT_INT_L2_MISS) | 
		                 EVNTSEL_SET(UNIT, EVNT_UNIT_TLB_RELOAD) |
		                 EVNTSEL_SET(USR, 1)                        |
		                 EVNTSEL_SET(EDGE_DETECT, 1),
		.desc          = "L2 internal misses (TLB)"
	}
};

static const int __evnts_selected[] = {
	EVENT_DC_REFILLS_SYS_ALL,
	EVENT_DC_REFILLS_L2_ALL,
        EVENT_L2_REQS_TLB_RELOAD,
	EVENT_L2_MISSES_TLB_RELOAD,
};

/*
 * those work only for xeon 
 */
#define PRFCNT_FL_T0 0 
#define PRFCNT_FL_T1 0

#include "prfcnt_simple.h"

#endif
