/*
 * spmv_loops.c
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdlib.h>
#include <sched.h>

#include "spmv_loops.h"
#include "vector.h"
#include "spmv_method.h"

#include "tsc.h"
#ifdef SPMV_PRFCNT
#   include "prfcnt.h"
#endif

float SPMV_NAME(_bench_loop) (SPMV_NAME(_fn_t) *fn, void *matrix,
                              unsigned long loops, unsigned long rows_nr,
                              unsigned long cols_nr)
{
	VECTOR_TYPE *x,*y;
	unsigned long i;
	tsc_t tsc;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);

	tsc_init(&tsc);
#ifdef SPMV_PRFCNT
	int cpu = 0;
	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);
	CPU_SET(cpu,&cpu_set);
	sched_setaffinity(getpid(), sizeof(cpu_set_t), &cpu_set);
	prfcnt_t prfcnt;
	prfcnt_init(&prfcnt, cpu,PRFCNT_FL_T0|PRFCNT_FL_T1);
#endif
	for (i = 0; i < loops; i++) {
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);
		VECTOR_NAME(_init)(y, (ELEM_TYPE) 0);
#ifdef SPMV_PRFCNT
		prfcnt_start(&prfcnt);
#endif
		tsc_start(&tsc);
		fn(matrix, x, y);
		tsc_pause(&tsc);
#ifdef SPMV_PRFCNT
		prfcnt_pause(&prfcnt);
#endif
	}

	tsc_shut(&tsc);
	const float secs = tsc_getsecs(&tsc);
#ifdef SPMV_PRFCNT
	prfcnt_report(&prfcnt);
	prfcnt_shut(&prfcnt);
#endif

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	return secs;
}

void SPMV_NAME(_check_loop) (void *m1, void *m2, SPMV_NAME(_fn_t) *fn1,
                             SPMV_NAME(_fn_t) *fn2, unsigned long loops,
                             unsigned long rows_nr, unsigned long cols_nr,
                             unsigned long nz_nr)
{
	unsigned long i;
	VECTOR_TYPE *x,*y1,*y2;

	x = VECTOR_NAME(_create)(cols_nr);
	y1 = VECTOR_NAME(_create)(rows_nr);
	y2 = VECTOR_NAME(_create)(rows_nr);

	for (i = 0; i < loops;  i++){
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);
		VECTOR_NAME(_init)(y1, (ELEM_TYPE) 23);
		VECTOR_NAME(_init)(y2, (ELEM_TYPE) 0);
		fn1(m1, x, y1);
		fn2(m2, x, y2);
		if (VECTOR_NAME(_compare)(y1, y2) < 0)
			exit(1);
	}

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y1);
	VECTOR_NAME(_destroy)(y2);
}
