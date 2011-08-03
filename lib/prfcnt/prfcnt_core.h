/*
 * prfcnt_core.h -- prfcnt implementation for Intel Core
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011,      Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __PRFCNT_CORE__
#define __PRFCNT_CORE__

#include "prfcnt_common.h"

/*
 * Pre-define Architectural (cross-platform) Performance Events
 * for Intel processors.
 *
 * For now We place them in the same file with the non-architectural
 * (Core - Specific) Events.
 *
 * NOTE: One can use the cpuid instruction to determine at
 * runtime things as number of counters per logical
 * processor, architectural events that are available
 * etc. Check cpuid.c.
 *
 * Here the definitions will be static.
 */

#define PRFCNT_EVNT_BADDR 0x186 /* event selectors base address */
#define PRFCNT_CNTR_BADDR 0x0c1 /* event counter base address */
#define PRFCNT_UNCORE_BADDR 0x107cc /* uncore event counter base address */

/*
 * selector registers fields
 * (duplicated from prfcnt_amd.h)
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

#define EVNT_UNIT_CORE_ALL       ((1<<6) | (1<<7))
#define EVNT_UNIT_CORE_THIS      ((1<<6) | (0<<7))
#define EVNT_UNIT_AGENT_ALL      (1<<5)
#define EVNT_UNIT_AGENT_THIS     (0<<5)
#define EVNT_UNIT_PFETCH_ALL     ((1<<4) | (1<<5))
#define EVNT_UNIT_PFETCH_HW      ((1<<4) | (0<<5))
#define EVNT_UNIT_PFETCH_NOHW    ((0<<4) | (0<<5))
#define EVNT_UNIT_MESI_INVALID   (1<<0)
#define EVNT_UNIT_MESI_SHARED    (1<<1)
#define EVNT_UNIT_MESI_EXCLUSIVE (1<<2)
#define EVNT_UNIT_MESI_MODIFIED  (1<<3)
#define EVNT_UNIT_MESI_ALL       ( EVNT_UNIT_MESI_INVALID   |\
                                   EVNT_UNIT_MESI_SHARED    |\
                                   EVNT_UNIT_MESI_EXCLUSIVE |\
                                   EVNT_UNIT_MESI_MODIFIED )

enum {
	EVNT_TYPE_PROGRAMMABLE, /* programmable event */
	EVNT_TYPE_FIXED,        /* fixed event */
	EVNT_TYPE_UNCORE,       /* event of the uncore of the
	                         * Intel Netburst/Core architectures */
};

struct prfcnt_core_event {
	__u64           evntsel_data;
	__u32           cntr_nr;
	unsigned long   flags;
	char            *desc;
	int             type;
};
typedef struct prfcnt_core_event prfcnt_event_t;

enum {
	PRFCNT_CORE_FL_ARCH_SHIFT = 0,
};
struct prfcnt_core_handle {
	int             msr_fd;
	unsigned int    cpu;
	unsigned long   flags;
	prfcnt_event_t  **events;
	unsigned int    events_nr;
};
typedef struct prfcnt_core_handle prfcnt_t;

/*
 * Available Events:
 * To add an Event you must first add a symbolic name
 * for it in the enum structure and then add the appropriate
 * details in __evnts[] array.
 * Events Reference:
 *  IA32 Software Developer's Manual: Volume 3B
 *  Performance Monitoring Events for INTEL CORE and INTEL CORE
 *  DUO PROCESSORS
 *  (Tables A-{8,9} )
 */

enum {
	EVENT_UNHLT_CORE_CYCLES = 0,
	EVENT_INSTR_RETIRED,
	EVENT_UOPS_RETIRED,
	EVENT_UOPS_ISSUED_ANY,
	EVENT_LLC_REFERENCE,
	EVENT_LLC_MISSES,
	EVENT_BR_RETIRED,
	EVENT_MISPRED_BR_RETIRED,
	EVENT_MUL,
	EVENT_DIV,
	EVENT_L1I_READS,
	EVENT_L1I_MISSES,
	EVENT_L1D_REPL,
	EVENT_L1D_ALL_REF,
	EVENT_L2_REJECT_CYCLES,
	EVENT_L2_LINES_IN,
	EVENT_L2_RQSTS,
	EVENT_RESOURCE_STALL,
	EVENT_FP_COMP_OPS_EXE,
	EVENT_BR_BOGUS,
	EVENT_BP_BTB_MISSES,
	EVENT_FUSED_UOPS_RET,
	EVENT_BR_IDEC,
	EVENT_BR_EXECUTED,
	EVENT_BR_MISSPREDICTED,
	EVENT_SIMD_UOPS_EXEC,
	EVENT_SIMD_UOP_TYPE_EXEC_MUL,
	EVENT_BUS_TRANS_MEM,
	EVENT_BUS_TRANS_ANY,
	EVENT_UNCORE_FSB_ALL,
	EVENT_END
};

static prfcnt_event_t __evnts[] = {
	/*
	 * Architectural Events
	 */
	[EVENT_UNHLT_CORE_CYCLES] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x3c)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Unhalted Core Cycles",
	},

	[EVENT_INSTR_RETIRED] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xc0)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Instructions Retired",
	},

	[EVENT_LLC_REFERENCE] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x2e)    |
		                EVNTSEL_SET(UNIT, 0x4f)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Last Level Cache references",
	},

	[EVENT_LLC_MISSES] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x2e)    |
		                EVNTSEL_SET(UNIT, 0x41)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Last Level cache misses",
	},

	[EVENT_BR_RETIRED] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xc4)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.desc         = "Branch Instruction Retired",
		.type         = EVNT_TYPE_PROGRAMMABLE,
	},

	[EVENT_MISPRED_BR_RETIRED] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xc5)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Mispredicted Branch Instruction Retired",
	},

	/*
	 * Non-Architectural Events
	 */
	[EVENT_UOPS_ISSUED_ANY] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x0e)    |
		                EVNTSEL_SET(UNIT, 0x01)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Uops issued",
	},

	[EVENT_UOPS_RETIRED] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xc2)    |
		                EVNTSEL_SET(UNIT, 0x0f)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Micro-ops Retired",
	},

	[EVENT_FP_COMP_OPS_EXE] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x10)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.desc         = "Floating point computational micro-ops executed",
		.type         = EVNT_TYPE_PROGRAMMABLE,
	},

	[EVENT_SIMD_UOPS_EXEC] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xb0)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "SIMD micro-ops executed",
	},

	[EVENT_SIMD_UOP_TYPE_EXEC_MUL] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xb3)    |
		                EVNTSEL_SET(UNIT, 0x01)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "SIMD packed multiply micro-ops executed",
	},

	[EVENT_MUL] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x12)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type          = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "Multiply Operations (INT+FP)",
	},

	[EVENT_DIV] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x13)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Divition Operations (INT+FP)",
	},

	[EVENT_L1I_READS] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x80)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Instruction fetches",
	},

	[EVENT_L1I_MISSES] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x81)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Instruction Fetch Unit misses",
	},

	[EVENT_L1D_ALL_REF] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x43)    |
		                EVNTSEL_SET(UNIT, 0x01)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "All references to the L1 data cache",
	},

	[EVENT_L1D_REPL] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x45)    |
		                EVNTSEL_SET(UNIT, 0x0f)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Cache lines allocated in the L1 data cache",
	},

	[EVENT_L2_REJECT_CYCLES] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x30)                |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_MESI_ALL)   |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_PFETCH_ALL) |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_CORE_THIS)  |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Cycles L2 is Busy (MESI_ALL,PF_ALL,CORE_THIS)",
	},

	[EVENT_L2_LINES_IN] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x24)                |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_CORE_THIS)  |
		                EVNTSEL_SET(USR, 1),
		.type          = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "L2 cache misses (no hwpref)",
	},

	[EVENT_L2_RQSTS] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x2e)                |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_MESI_ALL)   |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_CORE_THIS)  |
		                EVNTSEL_SET(USR, 1),
		.type          = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "L2 cache reference requests (no hwpref)",
	},

	[EVENT_RESOURCE_STALL] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xa2)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type          = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "Cycles of Resource Stall",
	},

	[EVENT_BR_BOGUS] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xe4)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type          = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "Bogus brances",
	},

	[EVENT_BP_BTB_MISSES] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xe2)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type          = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "branch misses in BTB",
	},

	[EVENT_BR_IDEC] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xe0)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type          = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "Branch instructions decoded",
	},

	[EVENT_FUSED_UOPS_RET] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0xda)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Fused uops retired"
	},

	[EVENT_BR_EXECUTED] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x88)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Branch instruction executed",
	},

	[EVENT_BR_MISSPREDICTED] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x89)    |
		                EVNTSEL_SET(UNIT, 0x00)     |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc          = "Branch instruction executed and mispredicted",
	},

	[EVENT_BUS_TRANS_MEM] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x6f)                |
		                EVNTSEL_SET(UNIT, 0xc0)                 |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_AGENT_ALL)  |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "Completed memory transactions",
	},

	[EVENT_BUS_TRANS_ANY] = {
		.evntsel_data = EVNTSEL_SET(EVENT, 0x70)                |
		                EVNTSEL_SET(UNIT, 0xc0)                 |
		                EVNTSEL_SET(UNIT, EVNT_UNIT_AGENT_ALL)  |
		                EVNTSEL_SET(USR, 1),
		.type         = EVNT_TYPE_PROGRAMMABLE,
		.desc         = "All completed bus transactions",
	},

	[EVENT_UNCORE_FSB_ALL] = {
		.cntr_nr      = 0,
		.evntsel_data = 0xff7L << 32,
		.type         = EVNT_TYPE_UNCORE,
		.desc         = "All bus transactions",
	},
};

/*
 * Select Events Here
 * (The number of available events can be determined using the cpuid
 * intruction)
 */
static const int __evnts_selected[] = {
/* 	EVENT_UNHLT_CORE_CYCLES, */
/* 	EVENT_UOPS_RETIRED, */
/* 	EVENT_L1D_REPL, */
/* 	EVENT_L1D_ALL_REF, */
/* 	EVENT_L1I_READS, */
/* 	EVENT_L1I_MISSES, */
/* 	EVENT_L2_LINES_IN, */
/* 	EVENT_L2_RQSTS, */
/* 	EVENT_LLC_REFERENCE, */
/* 	EVENT_LLC_MISSES, */
/* 	EVENT_BUS_TRANS_ANY, */
/* 	EVENT_SIMD_UOPS_EXEC, */
	EVENT_BR_RETIRED,
	EVENT_MISPRED_BR_RETIRED,
/*	EVENT_UNCORE_FSB_ALL, */
};

/*
 * those work only for xeon
 */
#define PRFCNT_FL_T0 0
#define PRFCNT_FL_T1 0

#include "prfcnt_simple.h"

#endif /* __PRFCNT_CORE__ */
