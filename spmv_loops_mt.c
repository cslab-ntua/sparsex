#include <stdlib.h>
#include <pthread.h>

#include "spmv_method.h"
#include "vector.h"
#include "spm_mt.h"
#include "mt_lib.h"

#include "spmv_loops_mt.h"
#include "tsc.h"

static VECTOR_TYPE *x=NULL, *y=NULL;
static pthread_barrier_t barrier;
static unsigned long loops_nr = 0;
static float secs = 0.0;

static void *do_spmv_thread(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *)arg;
	//printf("Hello I'm thread on cpu:%d\n", spm_mt_thread->cpu);
	SPMV_NAME(_fn_t) *spmv_mt_fn = spm_mt_thread->spmv_fn;
	setaffinity_oncpu(spm_mt_thread->cpu);
	int i;

	for (i=0; i<loops_nr; i++){
		pthread_barrier_wait(&barrier);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
		pthread_barrier_wait(&barrier);
	}

	return NULL;
}

#define SWAP(x,y) do {    \
	typeof(x) _tmp;   \
	_tmp = x;         \
	x = y;            \
	y = _tmp;         \
} while (0)

static void *do_spmv_thread_main(void *arg)
{
	spm_mt_thread_t *spm_mt_thread;
	SPMV_NAME(_fn_t) *spmv_mt_fn;
	tsc_t tsc;

	spm_mt_thread = arg;
	spmv_mt_fn = spm_mt_thread->spmv_fn;
	setaffinity_oncpu(spm_mt_thread->cpu);

	VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);

	//printf("Hello I'm thread on cpu:%d\n", spm_mt_thread->cpu);
	tsc_init(&tsc);
	tsc_start(&tsc);
	int i;
	for (i=0; i<loops_nr; i++){
		pthread_barrier_wait(&barrier);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
		SWAP(x, y);
		VECTOR_NAME(_init)(y, (ELEM_TYPE)0);
		pthread_barrier_wait(&barrier);
	}
	tsc_pause(&tsc);
	secs = tsc_getsecs(&tsc);
	tsc_shut(&tsc);
	return NULL;
}

static void init_barrier(unsigned count)
{
	int err;
	err = pthread_barrier_init(&barrier, NULL, count);
	if (err){
		perror("pthread_barrier_init");
		exit(1);
	}
}

float SPMV_NAME(_bench_mt_loop) (spm_mt_t *spm_mt, unsigned long loops,
				 unsigned long rows_nr, unsigned long cols_nr,
				 SPMV_NAME(_fn_t) *fn)
{
	int i;
	pthread_t *tids;

	secs = 0.0;
	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);
	loops_nr = loops;
	init_barrier(spm_mt->nr_threads);

	tids = malloc(sizeof(pthread_t)*spm_mt->nr_threads);
	if ( !tids ){
		perror("malloc");
		exit(1);
	}

	for (i=0; i<spm_mt->nr_threads; i++){
		if (fn != NULL)
			spm_mt->spm_threads[i].spmv_fn = fn;
	}

	pthread_create(tids, NULL, do_spmv_thread_main, spm_mt->spm_threads);
	for (i=1; i<spm_mt->nr_threads; i++)
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	for (i=0; i<spm_mt->nr_threads; i++)
		pthread_join(tids[i], NULL);

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	free(tids);

	return secs;
}

void SPMV_NAME(_check_mt_loop) ( void *spm, spm_mt_t *spm_mt,
                                 SPMV_NAME(_fn_t) *fn, unsigned long loops,
                                 unsigned long rows_nr, unsigned long cols_nr,
				 SPMV_NAME(_fn_t) *mt_fn)
{
	int err, i;
	pthread_t *tids;
	VECTOR_TYPE *y2;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);
	y2 = VECTOR_NAME(_create)(rows_nr);
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
		if (mt_fn != NULL)
			spm_mt->spm_threads[i].spmv_fn = mt_fn;
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);
	}

	for (i=0; i<loops; i++){
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);
		VECTOR_NAME(_init)(y, (ELEM_TYPE)0);
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

	for (i=0; i<spm_mt->nr_threads; i++){
		pthread_join(tids[i], NULL);
	}
	free(tids);
}

void SPMV_NAME(_check_mt_loop_serial) ( void *spm, spm_mt_t *spm_mt,
                                        SPMV_NAME(_fn_t) *fn, unsigned long loops,
                                        unsigned long rows_nr, unsigned long cols_nr,
				        SPMV_NAME(_fn_t) *mt_fn)
{
	int i, j;
	VECTOR_TYPE *y2;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);
	y2 = VECTOR_NAME(_create)(rows_nr);
	loops_nr = loops;

	for (i=0; i<spm_mt->nr_threads; i++){
		if (mt_fn != NULL)
			spm_mt->spm_threads[i].spmv_fn = mt_fn;
	}

	for (i=0; i<loops; i++){
		//unsigned long  k;
		//for (k=0; k < cols_nr; k++)
		//	x->elements[k] = (ELEM_TYPE)(k+666.0);
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);
		VECTOR_NAME(_init)(y, (ELEM_TYPE)0);
		VECTOR_NAME(_init)(y2, (ELEM_TYPE)21);

		for (j=0; j<spm_mt->nr_threads; j++){
			spm_mt_thread_t *t = spm_mt->spm_threads + j;
			SPMV_NAME(_fn_t) *spmv_mt_fn = t->spmv_fn;
			spmv_mt_fn(t->spm, x, y);
		}

		fn(spm, x, y2);

		if ( VECTOR_NAME(_compare)(y2, y) < 0){
			exit(1);
		}
	}

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	VECTOR_NAME(_destroy)(y2);
}
