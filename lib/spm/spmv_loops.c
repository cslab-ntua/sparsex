#include <stdlib.h>
#include <sched.h>

#include "spmv_loops.h"
#include "vector.h"
#include "spmv_method.h"

#include "tsc.h"
#ifdef SPMV_PRFCNT
#include "prfcnt.h"
#endif

//#define SPMV_PRFCNT

float SPMV_NAME(_bench_loop) (SPMV_NAME(_fn_t) *fn, void *matrix,
                             unsigned long loops,
			     unsigned long rows_nr, unsigned long cols_nr)
{

	VECTOR_TYPE *x,*y;
	unsigned long i;
	tsc_t tsc;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);

	tsc_init(&tsc);
	#ifdef SPMV_PRFCNT
	int cpu=0;
	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);
	CPU_SET(cpu,&cpu_set);
	sched_setaffinity(getpid(), sizeof(cpu_set_t), &cpu_set);
	prfcnt_t prfcnt;
	prfcnt_init(&prfcnt, cpu,PRFCNT_FL_T0|PRFCNT_FL_T1);
	#endif
	for (i=0; i<loops; i++){
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);
		VECTOR_NAME(_init)(y, (ELEM_TYPE)0);
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
	//const float mf = (((float)(nz_nr*2*loops))/secs)/(1000*1000);
	//const float mb = (((float)(nz_nr*2+rows_nr*2+cols_nr)*8*LOOPS)/secs)/(1024.0*1024.0);
	//printf("%s %f %f\n", prefix, mf, secs);
	#ifdef SPMV_PRFCNT
	prfcnt_report(&prfcnt);
	prfcnt_shut(&prfcnt);
	#endif

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);

	return secs;
}

void SPMV_NAME(_check_loop) (void *m1, void *m2,
                             SPMV_NAME(_fn_t) *fn1, SPMV_NAME(_fn_t) *fn2,
                             unsigned long loops,
                             unsigned long rows_nr, unsigned long cols_nr,
			     unsigned long nz_nr)
{
	unsigned long i;
	VECTOR_TYPE *x,*y1,*y2;

	x = VECTOR_NAME(_create)(cols_nr);
	y1 = VECTOR_NAME(_create)(rows_nr);
	y2 = VECTOR_NAME(_create)(rows_nr);

	for (i=0; i<loops ;  i++){
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);
		//vector_init(x, (ELEM_TYPE)3);
		//printf("y1=%p y2=%p\n", y1->elements, y2->elements);
		// We need to zero the vectors, for spmv versions that use
		// vertical or diagonal patterns
		VECTOR_NAME(_init)(y1, (ELEM_TYPE)0);
		VECTOR_NAME(_init)(y2, (ELEM_TYPE)0);
		fn1(m1, x, y1);
		fn2(m2, x, y2);
		if ( VECTOR_NAME(_compare)(y1, y2) < 0){
			exit(1);
		}
		//vector_print(x);
	}

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y1);
	VECTOR_NAME(_destroy)(y2);
}
