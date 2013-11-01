/*
 * CxSpmvMt.cc
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "CsxSpmvMt.hpp"

static vector_t *x = NULL;
static vector_t *y = NULL;
static vector_t **temp = NULL;
static unsigned long nr_loops = 0;
static unsigned long nr_threads = 0;
static unsigned long n = 0;
static float secs = 0.0;

static void do_spmv_thread(spm_mt_thread_t *thread, boost::barrier &cur_barrier)
{
	setaffinity_oncpu(thread->cpu);

	spmv_double_fn_t *spmv_mt_fn = (spmv_double_fn_t *) thread->spmv_fn;
    tsc_t thread_tsc;

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) thread->data;
	prfcnt_init(prfcnt, thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif

    tsc_init(&thread_tsc);
	for (unsigned int i = 0; i < nr_loops; i++) {
		cur_barrier.wait();
        tsc_start(&thread_tsc);
		spmv_mt_fn(thread->spm, x, y, 1);
        tsc_pause(&thread_tsc);
		cur_barrier.wait();
	}

#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif

    thread->secs = tsc_getsecs(&thread_tsc);
    tsc_shut(&thread_tsc);
}

#define SWAP(x,y) \
	do { \
		__typeof__(x) _tmp; \
		_tmp = x; \
		x = y; \
		y = _tmp; \
	} while (0)

static void do_spmv_thread_main_swap(spm_mt_thread_t *thread,
                                     boost::barrier &cur_barrier)
{
	setaffinity_oncpu(thread->cpu);

	spmv_double_fn_t *spmv_mt_fn = (spmv_double_fn_t *) thread->spmv_fn;
	tsc_t total_tsc, thread_tsc;

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) thread->data;
#endif

	vec_init_rand_range(x, (value_t) -1000, (value_t) 1000);

	// Assert this is a square matrix and swap is ok.
	assert(x->size == y->size);
	tsc_init(&total_tsc);
	tsc_init(&thread_tsc);
	tsc_start(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_init(prfcnt, thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif

	for (unsigned int i = 0; i < nr_loops; i++) {
		cur_barrier.wait();
		tsc_start(&thread_tsc);
		spmv_mt_fn(thread->spm, x, y, 1);
		tsc_pause(&thread_tsc);
		cur_barrier.wait();
		SWAP(x, y);
	}

	tsc_pause(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif
	thread->secs = tsc_getsecs(&thread_tsc);
	secs = tsc_getsecs(&total_tsc);
	tsc_shut(&total_tsc);
	tsc_shut(&thread_tsc);
}

static void do_spmv_thread_sym(spm_mt_thread_t *thread, boost::barrier &cur_barrier)
{
	spmv_double_sym_fn_t *spmv_mt_sym_fn = 
        (spmv_double_sym_fn_t *) thread->spmv_fn;
	int id = thread->id;
    tsc_t thread_tsc;
	// Switch Reduction Phase.
	/* int j, start, end*/;
	setaffinity_oncpu(thread->cpu);

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) thread->data;
	prfcnt_init(prfcnt, thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif
	// Switch Reduction Phase.
	/*
	start = (id * n) / ncpus;
	end = ((id + 1) * n) / ncpus;
	*/

    tsc_init(&thread_tsc);
	for (unsigned int i = 0; i < nr_loops; i++) {
		// Switch Reduction Phase.
		/*
		if (id != 0)
			vec_init(temp[id], 0);
			*/
		vec_init_from_map(temp, 0, thread->map);
		cur_barrier.wait();
        tsc_start(&thread_tsc);
		spmv_mt_sym_fn(thread->spm, x, y, temp[id]);
        tsc_pause(&thread_tsc);
		cur_barrier.wait();
		// Switch Reduction Phase.
		/*
		for (j = 1; j < ncpus; j++)
			vec_add_part(y, temp[j], y, start, end);
		*/
		vec_add_from_map(y, temp, y, thread->map);
		cur_barrier.wait();
	}

#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif
    thread->secs = tsc_getsecs(&thread_tsc);
    tsc_shut(&thread_tsc);
}

static void do_spmv_thread_sym_main_swap(spm_mt_thread_t *thread,
                                         boost::barrier &cur_barrier)
{
	setaffinity_oncpu(thread->cpu);
	spmv_double_sym_fn_t *spmv_mt_sym_fn = 
        (spmv_double_sym_fn_t *) thread->spmv_fn;
	tsc_t tsc, total_tsc, thread_tsc;
	// int id = thread->id;

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) thread->data;
#endif

	vec_init_rand_range(x, (value_t) -1000, (value_t) 1000);

	// Assert that the matrix is square and swap is OK.
	assert(x->size == y->size);
	tsc_init(&tsc);
	tsc_start(&tsc);
	tsc_init(&total_tsc);
	tsc_init(&thread_tsc);
	tsc_start(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_init(prfcnt, thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif
	// Switch Reduction Phase.
	/*int j, start, end*/;
	/*
	start = (id * n) / ncpus;
	end = ((id + 1) * n) / ncpus;
	*/
	for (unsigned int i = 0; i < nr_loops; i++) {
		// Switch Reduction Phase
		vec_init_from_map(temp, 0, thread->map);
		cur_barrier.wait();
		tsc_start(&thread_tsc);
		spmv_mt_sym_fn(thread->spm, x, y, y);
		tsc_pause(&thread_tsc);
		cur_barrier.wait();
		// Switch Reduction Phase
		/*
		for (j = 1; j < ncpus; j++)
			VECTOR_NAME(_add_part)(y, temp[j], y, start, end);
		*/
		vec_add_from_map(y, temp, y, thread->map);
		cur_barrier.wait();
		SWAP(x, y);
	}
	tsc_pause(&tsc);
	tsc_pause(&total_tsc);

#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif

	secs = tsc_getsecs(&tsc);
	tsc_shut(&tsc);
	thread->secs = tsc_getsecs(&thread_tsc);
	secs = tsc_getsecs(&total_tsc);
	tsc_shut(&total_tsc);
	tsc_shut(&thread_tsc);
}

static void csr_spmv(csx::CSR<index_t, value_t> *spm, vector_t *in, vector_t *out)
{
	value_t *x = in->elements;
	value_t *y = out->elements;
	value_t *values = spm->values_;
	index_t *row_ptr = spm->rowptr_;
	index_t *col_ind = spm->colind_;
	const index_t row_start = 0;
	const index_t row_end = spm->GetNrRows();
	register value_t yr;
	index_t i,j;

	for (i = row_start; i < row_end; i++) {
		yr = (value_t) 0;
		for (j = row_ptr[i]; j < row_ptr[i+1]; j++)
			yr += (values[j] * x[col_ind[j]]);

		y[i] = yr;
	}
}

// void csr_spmv_sym(csx::CSR<index_t, value_t> *spm, vector_t *in, vector_t *out)
// {
// 	value_t *y = out->elements;
// 	value_t *x = in->elements;
// 	value_t *values = spm->values_;
// 	value_t *dvalues = spm->dvalues_;
// 	index_t *row_ptr = spm->rowptr_;
// 	index_t *col_ind = spm->colind_;
// 	size_t n = spm->GetNrRows();
// 	register value_t yr;
// 	size_t i, j;

// 	for (i = 0; i < n; i++) {
// 		yr  = (value_t) 0;
// 		for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
// 			yr += values[j] * x[col_ind[j]];
// 			y[col_ind[j]] += values[j] * x[i];
// 		}

// 		yr += dvalues[i] * x[i];
// 		y[i] = yr;
// 	}
// }

float spmv_bench_mt(spm_mt_t *spm_mt, unsigned long loops,
                    unsigned long rows_nr, unsigned long cols_nr)
{
	secs = 0.0;
	nr_loops = loops;
    nr_threads = spm_mt->nr_threads;

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
		xparts[i] = spm->nr_rows * sizeof(value_t);
		yparts[i] = spm->nr_rows * sizeof(value_t);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
    }

	// Allocate an interleaved x.
	int alloc_err = 0;
	x = vec_create_interleaved(rows_nr, xparts, nr_threads, xnodes);
	vec_init_rand_range(x, (value_t) -1000, (value_t) 1000);
	alloc_err = check_interleaved(x->elements, xparts, nr_threads, xnodes);
	print_alloc_status("input vector", alloc_err);

	// Allocate an interleaved y.
	y = vec_create_interleaved(rows_nr, yparts, nr_threads, ynodes);
	vec_init(y, 0);
	alloc_err = check_interleaved(y->elements, yparts, nr_threads, ynodes);
	print_alloc_status("output vector", alloc_err);
#else
	x = vec_create(cols_nr, NULL);
	y = vec_create(rows_nr, NULL);
#endif

	/*
	 * spawn two kind of threads:
	 *	- 1	: do_spmv_thread_main_swap -> computes and does swap.
	 *	- N-1	: do_spmv_thread -> just computes.
	 */
    boost::barrier cur_barrier(nr_threads);
    vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
    threads[0] = boost::make_shared<boost::thread>
        (boost::bind(&do_spmv_thread_main_swap, spm_mt->spm_threads,
                         boost::ref(cur_barrier)));
    for (size_t i = 1; i < nr_threads; ++i) {
        threads[i] = boost::make_shared<boost::thread>
            (boost::bind(&do_spmv_thread,
                         spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

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
#ifdef SPM_NUMA
	free(xparts);
	free(yparts);
	free(xnodes);
	free(ynodes);
#endif

	return secs;
}

void spmv_check_mt(csx::CSR<index_t, value_t> *spm, spm_mt_t *spm_mt,
                   unsigned long loops, unsigned long rows_nr,
                   unsigned long cols_nr)
{
	vector_t *y2;
	nr_loops = loops;
    nr_threads = spm_mt->nr_threads;

#ifdef SPM_NUMA
	size_t *parts = (size_t *) malloc(sizeof(*parts)*spm_mt->nr_threads);
	int *nodes = (int *) malloc(sizeof(*nodes)*spm_mt->nr_threads);
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	for (size_t i = 0; i < spm_mt->nr_threads; i++) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
		parts[i] = spm->nr_rows * sizeof(value_t);
		nodes[i] = spm->node;
	}

	int node = spm_mt->spm_threads[0].node;
	x = vec_create_onnode(cols_nr, node);
	y = vec_create_interleaved(rows_nr, parts, spm_mt->nr_threads, nodes);
	y2 = vec_create(rows_nr, NULL);
#else
	x = vec_create(cols_nr, NULL);
	y = vec_create(rows_nr, NULL);
	y2 = vec_create(rows_nr, NULL);
#endif
	vec_init_rand_range(x, (value_t) -1000, (value_t) 1000);
	vec_init(y, (value_t) 0);
    vec_init(y2, (value_t) 0);

    vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
    boost::barrier cur_barrier(nr_threads + 1);

	for (unsigned int i = 0; i < nr_threads; ++i) {
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
        threads[i] = boost::make_shared<boost::thread>
            (boost::bind(&do_spmv_thread,
                         spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
	}

	for (unsigned int i = 0; i < loops; i++) {
		cur_barrier.wait();
		cur_barrier.wait();
		csr_spmv(spm, x, y2);
		if (vec_compare(y2, y) < 0)
			exit(1);
	}

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

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
#ifdef SPM_NUMA
	free(nodes);
	free(parts);
#endif
}

float spmv_bench_sym_mt(spm_mt_t *spm_mt, unsigned long loops,
                        unsigned long nrows, unsigned long ncols)
{
	secs = 0.0;
	nr_loops = loops;
	nr_threads = spm_mt->nr_threads;
	n = nrows;

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
		xparts[i] = spm->nr_rows * sizeof(value_t);
		yparts[i] = spm->nr_rows * sizeof(value_t);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
    }

	// Allocate an interleaved x.
	int alloc_err = 0;
	x = vec_create_interleaved(nrows, xparts, nr_threads, xnodes);
	vec_init_rand_range(x, (value_t) -1000, (value_t) 1000);
	alloc_err = check_interleaved(x->elements, xparts, nr_threads, xnodes);
	print_alloc_status("input vector", alloc_err);

	// Allocate an interleaved y.
	y = vec_create_interleaved(nrows, yparts, nr_threads, ynodes);
	vec_init(y, 0);
	alloc_err = check_interleaved(y->elements, yparts, nr_threads, ynodes);
	print_alloc_status("output vector", alloc_err);
#else
	x = vec_create(n, NULL);
	y = vec_create(n, NULL);
#endif

	temp = (vector_t **) malloc(nr_threads * sizeof(vector_t *));
	if (!temp) {
		perror("malloc");
		exit(1);
	}

    temp[0] = y;
	for (unsigned int i = 1; i < nr_threads; i++) {
#ifdef SPM_NUMA
		int tnode = spm_mt->spm_threads[i].node;
		temp[i] = vec_create_onnode(n, tnode);
#else
		temp[i] = vec_create(n, NULL);
#endif
    }

	for (unsigned int i = 1; i < nr_threads; i++) {
		for (unsigned int j = 1; j < n; j++)
			temp[i]->elements[j] = 0;
	}

#ifdef SPM_NUMA
	alloc_err = 0;
	for (unsigned int i = 1; i < nr_threads; i++) {
		int tnode = spm_mt->spm_threads[i].node;
		alloc_err += check_region(temp[i]->elements, n*sizeof(value_t), tnode);
	}
	print_alloc_status("local buffers", alloc_err);
#endif

    boost::barrier cur_barrier(nr_threads);
    vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
    threads[0] = boost::make_shared<boost::thread>
        (boost::bind(&do_spmv_thread_sym_main_swap, spm_mt->spm_threads,
                         boost::ref(cur_barrier)));
    for (size_t i = 1; i < nr_threads; ++i) {
        threads[i] = boost::make_shared<boost::thread>
            (boost::bind(&do_spmv_thread_sym,
                         spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

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
	for (unsigned int i = 1; i < nr_threads; i++)
		vec_destroy(temp[i]);
	free(temp);
#ifdef SPM_NUMA
    free(xparts);
    free(xnodes);
    free(yparts);
    free(ynodes);
#endif

	return secs;
}

void spmv_check_sym_mt(csx::CSR<index_t, value_t> *spm, spm_mt_t *spm_mt,
                       unsigned long loops, unsigned long nrows,
                       unsigned long ncols)
{
	vector_t *y2;

	nr_loops = loops;
	nr_threads = spm_mt->nr_threads;
	n = nrows;
//    assert(nrows == ncols);

#ifdef SPM_NUMA
	size_t *parts = (size_t *) malloc(nr_threads * sizeof(*parts));
	int *nodes = (int *) malloc(nr_threads * sizeof(*nodes));
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
		parts[i] = spm->nr_rows * sizeof(value_t);
		nodes[i] = spm->node;
	}

	x = vec_create_interleaved(n, parts, nr_threads, nodes);
	y = vec_create_interleaved(n, parts, nr_threads, nodes);
	y2 = vec_create(n, NULL);
#else
	x = vec_create(n, NULL);
	y = vec_create(n, NULL);
	y2 = vec_create(n, NULL);
#endif
    vec_init_rand_range(x, (value_t) -1000, (value_t) 1000);
    vec_init(y, (value_t) 0);
    vec_init(y2, (value_t) 0);

	temp = (vector_t **) malloc(nr_threads * sizeof(vector_t *));
	if (!temp) {
		perror("malloc");
		exit(1);
	}

	temp[0] = y;
	for (unsigned int i = 1; i < nr_threads; i++) {
#ifdef SPM_NUMA
		int node = spm_mt->spm_threads[i].node;
        temp[i] = vec_create_onnode(n, node);
		for (size_t j = 1; j < n; j++)
			temp[i]->elements[j] = 0;
#else
		temp[i] = vec_create(n, NULL);
#endif
    }

    vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
    boost::barrier cur_barrier(nr_threads + 1);

	for (unsigned int i = 0; i < nr_threads; i++) {
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
        threads[i] = boost::make_shared<boost::thread>
            (boost::bind(&do_spmv_thread_sym,
                         spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
	}

	for (unsigned int i = 0; i < loops; i++) {
		cur_barrier.wait();
		cur_barrier.wait();
		cur_barrier.wait();
		csr_spmv(spm, x, y2);
		if (vec_compare(y2, y) < 0)
			exit(1);
	}

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

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
	for (unsigned int i = 1; i < nr_threads; i++)
		vec_destroy(temp[i]);
	free(temp);
#ifdef SPM_NUMA
    free(parts);
    free(nodes);
#endif
}
