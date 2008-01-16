#include <stdlib.h>
#include <pthread.h>

#include "spmv_method.h"
#include "vector.h"
#include "spm_mt.h"
#include "mt_lib.h"

#include "tsc.h"

static SPMV_NAME(_fn_t) *spmv_mt_fn=NULL;
static VECTOR_TYPE *x=NULL, *y=NULL;
static pthread_barrier_t barrier;
static unsigned long loops_nr = 0;

static void *do_spmv_thread(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *)arg;
	//printf("Hello I'm thread on cpu:%d\n", spm_mt_thread->cpu);
	setaffinity_oncpu(spm_mt_thread->cpu);
	int i;

	for (i=0; i<loops_nr; i++){
		pthread_barrier_wait(&barrier);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
		pthread_barrier_wait(&barrier);
	}

	return NULL;
}

float SPMV_NAME(_bench_mt_loop) ( SPMV_NAME(_fn_t) *fn, spm_mt_t *spm_mt,
                                  unsigned long loops, unsigned long cols_nr)
{
	tsc_t tsc;
	int err, i;
	pthread_t *tids;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(cols_nr);
	spmv_mt_fn = fn;
	loops_nr = loops;

	err = pthread_barrier_init(&barrier, NULL, spm_mt->nr_threads + 1);
	if (err){
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = malloc(sizeof(pthread_t)*spm_mt->nr_threads);
	if ( !tids ){
		perror("malloc");
		exit(1);
	}

	for (i=0; i<spm_mt->nr_threads; i++){
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);
	}

	tsc_init(&tsc);
	for (i=0; i<loops; i++){
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);
		VECTOR_NAME(_init)(y, (ELEM_TYPE)0);

		tsc_start(&tsc);
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		tsc_pause(&tsc);
	}
	tsc_shut(&tsc);

	const float secs = tsc_getsecs(&tsc);

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	free(tids);

	return secs;
}

void SPMV_NAME(_check_mt_loop) ( void *spm, spm_mt_t *spm_mt,
                                 SPMV_NAME(_fn_t) *fn, SPMV_NAME(_fn_t) *mt_fn,
				 unsigned long loops, unsigned long cols_nr)
{
	int err, i;
	pthread_t *tids;
	VECTOR_TYPE *y2;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(cols_nr);
	y2 = VECTOR_NAME(_create)(cols_nr);
	spmv_mt_fn = mt_fn;
	loops_nr = loops;

	err = pthread_barrier_init(&barrier, NULL, spm_mt->nr_threads + 1);
	if (err){
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = malloc(sizeof(pthread_t)*spm_mt->nr_threads);
	if ( !tids ){
		perror("malloc");
		exit(1);
	}

	for (i=0; i<spm_mt->nr_threads; i++){
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);
	}

	for (i=0; i<loops; i++){
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);
		VECTOR_NAME(_init)(y, (ELEM_TYPE)13);
		VECTOR_NAME(_init)(y2, (ELEM_TYPE)21);
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		fn(spm, x, y2);

		if ( VECTOR_NAME(_compare)(y2, y) < 0){
			exit(1);
		}
	}

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	VECTOR_NAME(_destroy)(y2);
}
