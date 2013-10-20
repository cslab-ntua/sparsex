/*
 * csx_spmv_mt.cc
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "csx_spmv_mt.h"

static vector_t *x = NULL;
static vector_t *y = NULL;
static pthread_barrier_t barrier;
static unsigned long loops_nr = 0;
static float secs = 0.0;

static void *do_spmv_thread(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *) arg;
	spmv_double_fn_t *spmv_mt_fn = (spmv_double_fn_t *) spm_mt_thread->spmv_fn;
	setaffinity_oncpu(spm_mt_thread->cpu);
    tsc_t thread_tsc;

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) spm_mt_thread->data;
	prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif

    tsc_init(&thread_tsc);
	for (unsigned int i = 0; i < loops_nr; i++) {
		pthread_barrier_wait(&barrier);
        tsc_start(&thread_tsc);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
        tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
	}

#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif

    spm_mt_thread->secs = tsc_getsecs(&thread_tsc);
    tsc_shut(&thread_tsc);
	return NULL;
}

#define SWAP(x,y) \
	do { \
		__typeof__(x) _tmp; \
		_tmp = x; \
		x = y; \
		y = _tmp; \
	} while (0)

static void *do_spmv_thread_main_swap(void *arg)
{
	spm_mt_thread_t *spm_mt_thread;
#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt;
#endif
	spmv_double_fn_t *spmv_mt_fn;
	tsc_t total_tsc, thread_tsc;

	spm_mt_thread = (spm_mt_thread_t *) arg;
	spmv_mt_fn = (spmv_double_fn_t *) spm_mt_thread->spmv_fn;
#ifdef SPMV_PRFCNT
	prfcnt = (prfcnt_t *) spm_mt_thread->data;
#endif
	setaffinity_oncpu(spm_mt_thread->cpu);

	vec_init_rand_range(x, (double) -1000, (double) 1000);

	// Assert this is a square matrix and swap is ok.
	assert(x->size == y->size);
	tsc_init(&total_tsc);
	tsc_init(&thread_tsc);
	tsc_start(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif
	for (unsigned int i = 0; i < loops_nr; i++) {
		pthread_barrier_wait(&barrier);
		tsc_start(&thread_tsc);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
		tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
		SWAP(x, y);
	}
	tsc_pause(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif
	spm_mt_thread->secs = tsc_getsecs(&thread_tsc);
	secs = tsc_getsecs(&total_tsc);
	tsc_shut(&total_tsc);
	tsc_shut(&thread_tsc);

	return NULL;
}

float spmv_bench_mt(spm_mt_t *spm_mt, unsigned long loops,
                    unsigned long rows_nr, unsigned long cols_nr)
{
	secs = 0.0;
	loops_nr = loops;
    unsigned int nr_threads = spm_mt->nr_threads;

	int err = pthread_barrier_init(&barrier, NULL, nr_threads);
	if (err) {
		perror("pthread_barrier_init");
		exit(1);
	}

	pthread_t *tids = (pthread_t *) malloc(sizeof(pthread_t) * nr_threads);
	if (!tids) {
		perror("malloc");
		exit(1);
	}

#ifdef SPMV_PRFCNT
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
	}
#endif

#ifdef SPM_NUMA
	size_t *xparts = (size_t *) malloc(sizeof(*xparts) * nr_threads);
	size_t *yparts = (size_t *) malloc(sizeof(*yparts) * nr_threads);
	int *xnodes = (int *) malloc(sizeof(*xnodes) * nr_threads);
	int *ynodes = (int *) malloc(sizeof(*ynodes) * nr_threads);
	if (!xparts || !yparts || !xnodes || !ynodes) {
		perror("malloc");
		exit(1);
	}

	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm = spm_mt->spm_threads + i;
		xparts[i] = spm->nr_rows * sizeof(double);
		yparts[i] = spm->nr_rows * sizeof(double);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
    }

	// Allocate an interleaved x.
	int alloc_err = 0;
	x = vec_create_interleaved(rows_nr, xparts, nr_threads, xnodes);
	vec_init_rand_range(x, (double) -1000, (double) 1000);
	alloc_err = check_interleaved(x->elements, xparts, nr_threads, xnodes);
	print_alloc_status("input vector", alloc_err);

	// Allocate an interleaved y.
	y = vec_create_interleaved(rows_nr, yparts, nr_threads, ynodes);
	vec_init(y, 0);
	alloc_err = check_interleaved(y->elements, yparts, nr_threads, ynodes);
	print_alloc_status("output vector", alloc_err);
#else
	x = vec_create(cols_nr);
	y = vec_create(rows_nr);
#endif

	/*
	 * spawn two kind of threads:
	 *	- 1	: do_spmv_thread_main_swap -> computes and does swap.
	 *	- N-1	: do_spmv_thread -> just computes.
	 */
	pthread_create(tids, NULL, do_spmv_thread_main_swap, spm_mt->spm_threads);
	for (unsigned int i = 1; i < nr_threads; i++)
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	for (unsigned int i = 0; i < nr_threads; i++)
		pthread_join(tids[i], NULL);

#ifdef SPMV_PRFCNT
	// Report performance counters for every thread.
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
		prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
		fprintf(stdout, "Perf. counters: thread %d on cpu %d\n", i,
		        spmv_thread.cpu);
		prfcnt_report(prfcnt);
		prfcnt_shut(prfcnt);
		free(prfcnt);
	}
#endif

	vec_destroy(x);
	vec_destroy(y);
	pthread_barrier_destroy(&barrier);
#ifdef SPM_NUMA
	free(xparts);
	free(yparts);
	free(xnodes);
	free(ynodes);
#endif
	free(tids);
	return secs;
}

static void csr_multiply(csx::CSR<uint64_t, double> *spm, vector_t *in, vector_t *out)
{
	double *x = in->elements;
	double *y = out->elements;
	double *values = spm->values_;
	uint64_t *row_ptr = spm->rowptr_;
	uint64_t *col_ind = spm->colind_;
	const unsigned long row_start = 0;
	const unsigned long row_end = spm->GetNrRows();
	register double yr;
	unsigned long i,j;

	for (i = row_start; i < row_end; i++) {
		yr = (double) 0;
		for (j = row_ptr[i]; j < row_ptr[i+1]; j++)
			yr += (values[j] * x[col_ind[j]]);

		y[i] = yr;
	}
}

void spmv_check_mt(csx::CSR<uint64_t, double> *spm, spm_mt_t *spm_mt,
                   unsigned long loops, unsigned long rows_nr,
                   unsigned long cols_nr)
{
    unsigned int nr_threads = spm_mt->nr_threads;
	vector_t *y2;

	x = vec_create(cols_nr);
	y = vec_create(rows_nr);
	y2 = vec_create(rows_nr);
	loops_nr = loops;

	int err = pthread_barrier_init(&barrier, NULL, nr_threads + 1);
	if (err) {
		perror("pthread_barrier_init");
		exit(1);
	}

	pthread_t *tids = (pthread_t *) malloc(sizeof(pthread_t) * nr_threads);
	if (!tids) {
		perror("malloc");
		exit(1);
	}

	for (unsigned int i = 0; i < nr_threads; i++) {
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
		pthread_create(tids + i, NULL, do_spmv_thread, spm_mt->spm_threads + i);
	}

	for (unsigned int i = 0; i < loops; i++) {
		vec_init_rand_range(x, (double) -1000, (double) 1000);
		vec_init(y2, (double) 21);
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		csr_multiply(spm, x, y2);
		if (vec_compare(y2, y) < 0)
			exit(1);
	}

	for (unsigned int i = 0; i < nr_threads; i++)
		pthread_join(tids[i], NULL);

#ifdef SPMV_PRFCNT
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
		prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
		prfcnt_shut(prfcnt);
		free(prfcnt);
	}
#endif

	vec_destroy(x);
	vec_destroy(y);
	vec_destroy(y2);
	pthread_barrier_destroy(&barrier);
	free(tids);
}
