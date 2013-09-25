/*
 * csxsym_spmv_mt.cc
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "csxsym_spmv_mt.h"

static vector_t *x = NULL;
static vector_t *y = NULL;
static vector_t **temp = NULL;
static pthread_barrier_t barrier;
static unsigned long nloops = 0;
static unsigned long ncpus = 0;
static unsigned long n = 0;
static float secs = 0.0;

static void *do_spmv_thread(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *) arg;
	spmv_double_sym_fn_t *spmv_mt_sym_fn =
        (spmv_double_sym_fn_t *) spm_mt_thread->spmv_fn;
	int id = spm_mt_thread->id;
	// Switch Reduction Phase.
	/* int j, start, end*/;
    tsc_t thread_tsc;

	setaffinity_oncpu(spm_mt_thread->cpu);

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) spm_mt_thread->data;
	prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif
	// Switch Reduction Phase.
	/*
	start = (id * n) / ncpus;
	end = ((id + 1) * n) / ncpus;
	*/

    tsc_init(&thread_tsc);
	for (unsigned int i = 0; i < nloops; i++) {
		// Switch Reduction Phase.
		/*
		if (id != 0)
			vec_init(temp[id], 0);
			*/
		vec_init_from_map(temp, 0, spm_mt_thread->map);
		pthread_barrier_wait(&barrier);
        tsc_start(&thread_tsc);
		spmv_mt_sym_fn(spm_mt_thread->spm, x, y, temp[id]);
        tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
		// Switch Reduction Phase.
		/*
		for (j = 1; j < ncpus; j++)
			vec_add_part(y, temp[j], y, start, end);
		*/
		vec_add_from_map(y, temp, y, spm_mt_thread->map);
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
		typeof(x) _tmp;	\
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

	spmv_double_sym_fn_t *spmv_mt_sym_fn;
	tsc_t tsc;
	tsc_t total_tsc, thread_tsc;

	spm_mt_thread = (spm_mt_thread_t *) arg;
	spmv_mt_sym_fn = (spmv_double_sym_fn_t *) spm_mt_thread->spmv_fn;
	
	// int id = spm_mt_thread->id;

#ifdef SPMV_PRFCNT
	prfcnt = (prfcnt_t *) spm_mt_thread->data;
#endif
	setaffinity_oncpu(spm_mt_thread->cpu);

	vec_init_rand_range(x, (double) -1000, (double) 1000);

	// Assert that the matrix is square and swap is OK.
	assert(x->size == y->size);
	tsc_init(&tsc);
	tsc_start(&tsc);
	tsc_init(&total_tsc);
	tsc_init(&thread_tsc);
	tsc_start(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif
	// Switch Reduction Phase.
	/*int j, start, end*/;
	/*
	start = (id * n) / ncpus;
	end = ((id + 1) * n) / ncpus;
	*/
	for (unsigned int i = 0; i < nloops; i++) {
		// Switch Reduction Phase
		vec_init_from_map(temp, 0, spm_mt_thread->map);
		pthread_barrier_wait(&barrier);
		tsc_start(&thread_tsc);
		spmv_mt_sym_fn(spm_mt_thread->spm, x, y, y);
		tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
		// Switch Reduction Phase
		/*
		for (j = 1; j < ncpus; j++)
			VECTOR_NAME(_add_part)(y, temp[j], y, start, end);
		*/
		vec_add_from_map(y, temp, y, spm_mt_thread->map);
		pthread_barrier_wait(&barrier);
		SWAP(x, y);
	}
	tsc_pause(&tsc);
	tsc_pause(&total_tsc);

#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif

	secs = tsc_getsecs(&tsc);
	tsc_shut(&tsc);
	spm_mt_thread->secs = tsc_getsecs(&thread_tsc);
	secs = tsc_getsecs(&total_tsc);
	tsc_shut(&total_tsc);
	tsc_shut(&thread_tsc);

	return NULL;
}

float spmv_bench_sym_mt(spm_mt_t *spm_mt, unsigned long loops,
                        unsigned long nrows, unsigned long ncols)
{
	secs = 0.0;
	nloops = loops;
	ncpus = spm_mt->nr_threads;
	n = nrows;

	temp = (vector_t **) malloc(ncpus * sizeof(vector_t *));
	if (!temp) {
		perror("malloc");
		exit(1);
	}

	int err = pthread_barrier_init(&barrier, NULL, ncpus);
	if (err) {
		perror("pthread_barrier_init");
		exit(1);
	}

	pthread_t *tids = (pthread_t *) malloc(ncpus * sizeof(pthread_t));
	if (!tids) {
		perror("malloc");
		exit(1);
	}

#ifdef SPMV_PRFCNT
	for (unsigned int i = 0; i < ncpus; i++) {
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
	}
#endif

#ifdef SPM_NUMA
	size_t *xparts = (size_t *) malloc(sizeof(*xparts) * ncpus);
	size_t *yparts = (size_t *) malloc(sizeof(*yparts) * ncpus);
	int *xnodes = (int *) malloc(sizeof(*xnodes) * ncpus);
	int *ynodes = (int *) malloc(sizeof(*ynodes) * ncpus);
	if (!xparts || !yparts || !xnodes || !ynodes) {
		perror("malloc");
		exit(1);
	}

	for (unsigned int i = 0; i < ncpus; i++) {
		spm_mt_thread_t *spm = spm_mt->spm_threads + i;
		xparts[i] = spm->nr_rows * sizeof(double);
		yparts[i] = spm->nr_rows * sizeof(double);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
    }

	// Allocate an interleaved x.
	int alloc_err = 0;
	x = vec_create_interleaved(nrows, xparts, ncpus, xnodes);
	vec_init_rand_range(x, (double) -1000, (double) 1000);
	alloc_err = check_interleaved(x->elements, xparts, ncpus, xnodes);
	print_alloc_status("input vector", alloc_err);

	// Allocate an interleaved y.
	y = vec_create_interleaved(nrows, yparts, ncpus, ynodes);
	vec_init(y, 0);
	alloc_err = check_interleaved(y->elements, yparts, ncpus, ynodes);
	print_alloc_status("output vector", alloc_err);
#else
	x = vec_create(n);
	y = vec_create(n);
#endif
    temp[0] = y;
	for (unsigned int i = 1; i < ncpus; i++) {
#ifdef SPM_NUMA
		int tnode = spm_mt->spm_threads[i].node;
		temp[i] = vec_create_onnode(n, tnode);
#else
		temp[i] = vec_create(n);
#endif
    }

	for (unsigned int i = 1; i < ncpus; i++) {
		for (unsigned int j = 1; j < n; j++)
			temp[i]->elements[j] = 0;
	}

#ifdef SPM_NUMA
	alloc_err = 0;
	for (unsigned int i = 1; i < ncpus; i++) {
		int tnode = spm_mt->spm_threads[i].node;
		alloc_err += check_region(temp[i]->elements, n*sizeof(double),
		                          tnode);
	}
	print_alloc_status("local buffers", alloc_err);
#endif

	pthread_create(tids, NULL, do_spmv_thread_main_swap, spm_mt->spm_threads);
	for (unsigned int i = 1; i < ncpus; i++)
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	for (unsigned int i = 0; i < ncpus; i++)
		pthread_join(tids[i], NULL);

#ifdef SPMV_PRFCNT
	// Report performance counters for every thread.
	for (unsigned int i = 0; i < ncpus; i++) {
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
	for (unsigned int i = 1; i < ncpus; i++)
		vec_destroy(temp[i]);
	free(temp);
#ifdef SPM_NUMA
    free(xparts);
    free(xnodes);
    free(yparts);
    free(ynodes);
#endif
	pthread_barrier_destroy(&barrier);
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

void spmv_check_sym_mt(csx::CSR<uint64_t, double> *spm, spm_mt_t *spm_mt,
                       unsigned long loops, unsigned long nrows,
                       unsigned long ncols)
{
	vector_t *y2;

	nloops = loops;
	ncpus = spm_mt->nr_threads;
	n = nrows;

	x = vec_create(n);
	y = vec_create(n);
	y2 = vec_create(n);
	temp = (vector_t **) malloc(ncpus * sizeof(vector_t *));
	if (!temp) {
		perror("malloc");
		exit(1);
	}

	for (unsigned int i = 1; i < ncpus; i++)
		temp[i] = vec_create(n);
	temp[0] = y;
	
	int err = pthread_barrier_init(&barrier, NULL, ncpus+1);
	if (err) {
		perror("pthread_barrier_init");
		exit(1);
	}

	pthread_t *tids = (pthread_t *) malloc(ncpus * sizeof(pthread_t));
	if (!tids) {
		perror("malloc");
		exit(1);
	}

	for (unsigned int i = 0; i < ncpus; i++) {
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);
	}

	for (unsigned int i = 0; i < nloops; i++) {
		vec_init_rand_range(x, (double) -1000, (double) 1000);
		vec_init(y, (double) 21);
		vec_init(y2, (double) 22);
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		csr_multiply(spm, x, y2);
		if (vec_compare(y2, y) < 0)
			exit(1);
	}

	for (unsigned int i = 0; i < ncpus; i++)
		pthread_join(tids[i], NULL);

#ifdef SPMV_PRFCNT
	for (unsigned int i = 0; i < ncpus; i++) {
		spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
		prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
		prfcnt_shut(prfcnt);
		free(prfcnt);
	}
#endif

	vec_destroy(x);
	vec_destroy(y);
	vec_destroy(y2);
	for (unsigned int i = 1; i < ncpus; i++)
		vec_destroy(temp[i]);

	free(temp);
	pthread_barrier_destroy(&barrier);
	free(tids);
}
