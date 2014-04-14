/*
 * CsxBench.cpp -- Testing and benchmark utilities.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "CsxBench.hpp"

static spx_vector_t *x = NULL;
static spx_vector_t *y = NULL;
static spx_vector_t **tmp = NULL;
static size_t nr_loops = 0;
static size_t nr_threads = 0;
static size_t n = 0;
static float secs = 0.0;

#define SWAP(x,y) \
	do { \
		__typeof__(x) _tmp; \
		_tmp = x; \
		x = y; \
		y = _tmp; \
	} while (0)

int GetOptionOuterLoops()
{
    const char *loops_env = getenv("OUTER_LOOPS");
    int ret = 1;
    
    if (loops_env) {
        ret = atoi(loops_env);
        if (ret < 0)
            ret = 0;
    }
    
    return ret;
}

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
	for (size_t i = 0; i < nr_loops; i++) {
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

#ifdef DISABLE_POOL
static void do_spmv_thread_main_swap(spm_mt_thread_t *thread,
                                     boost::barrier &cur_barrier)
{
	setaffinity_oncpu(thread->cpu);
	spmv_double_fn_t *spmv_mt_fn = (spmv_double_fn_t *) thread->spmv_fn;
	tsc_t total_tsc, thread_tsc;

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) thread->data;
#endif

	vec_init_rand_range(x, (spx_value_t) -1000, (spx_value_t) 1000);

	// Assert this is a square matrix and swap is ok.
	assert(x->size == y->size);
	tsc_init(&total_tsc);
	tsc_init(&thread_tsc);
	tsc_start(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_init(prfcnt, thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif

	for (size_t i = 0; i < nr_loops; i++) {
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
#endif

static void do_spmv_thread_sym(spm_mt_thread_t *thread,
                               boost::barrier &cur_barrier)
{
	setaffinity_oncpu(thread->cpu);
	spmv_double_sym_fn_t *spmv_mt_sym_fn = 
        (spmv_double_sym_fn_t *) thread->spmv_fn;
	int id = thread->id;
    tsc_t thread_tsc;
	// Switch Reduction Phase.
	/* int j, start, end*/;

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
	for (size_t i = 0; i < nr_loops; i++) {
		// Switch Reduction Phase.
		/*
		if (id != 0)
			vec_init(tmp[id], 0);
			*/
		vec_init_from_map(tmp, 0, thread->map);
		cur_barrier.wait();
        tsc_start(&thread_tsc);
		spmv_mt_sym_fn(thread->spm, x, y, tmp[id], 1);
        tsc_pause(&thread_tsc);
		cur_barrier.wait();
		// Switch Reduction Phase.
		/*
		for (j = 1; j < ncpus; j++)
			vec_add_part(y, tmp[j], y, start, end);
		*/
		vec_add_from_map(y, tmp, y, thread->map);
		cur_barrier.wait();
	}

#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif
    thread->secs = tsc_getsecs(&thread_tsc);
    tsc_shut(&thread_tsc);
}

#ifdef DISABLE_POOL
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

	vec_init_rand_range(x, (spx_value_t) -1000, (spx_value_t) 1000);

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
	for (size_t i = 0; i < nr_loops; i++) {
		// Switch Reduction Phase
		vec_init_from_map(tmp, 0, thread->map);
		cur_barrier.wait();
		tsc_start(&thread_tsc);
		spmv_mt_sym_fn(thread->spm, x, y, y, 1);
		tsc_pause(&thread_tsc);
		cur_barrier.wait();
		// Switch Reduction Phase
		/*
		for (j = 1; j < ncpus; j++)
			VECTOR_NAME(_add_part)(y, tmp[j], y, start, end);
		*/
		vec_add_from_map(y, tmp, y, thread->map);
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
#endif

static void csr_spmv(csx::CSR<uindex_t, spx_value_t> *spm, 
                     spx_vector_t *in, spx_vector_t *out)
{
	spx_value_t *x = in->elements;
	spx_value_t *y = out->elements;
	spx_value_t *values = spm->values_;
	uindex_t *row_ptr = spm->rowptr_;
	uindex_t *col_ind = spm->colind_;
	const uindex_t row_start = 0;
	const uindex_t row_end = spm->GetNrRows();
	register spx_value_t yr;
	uindex_t i,j;

	for (i = row_start; i < row_end; i++) {
		yr = (spx_value_t) 0;
		for (j = row_ptr[i]; j < row_ptr[i+1]; j++)
			yr += (values[j] * x[col_ind[j]]);

		y[i] = yr;
	}
}

float spmv_bench_mt(spm_mt_t *spm_mt, size_t loops, size_t nr_rows,
                    size_t nr_cols)
{
	secs = 0.0;
	nr_loops = loops;
    nr_threads = spm_mt->nr_threads;

#ifdef SPMV_PRFCNT
	for (size_t i = 0; i < nr_threads; i++) {
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

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm = spm_mt->spm_threads + i;
		xparts[i] = spm->nr_rows * sizeof(spx_value_t);
		yparts[i] = spm->nr_rows * sizeof(spx_value_t);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
    }

	// Allocate an interleaved x.
	x = vec_create_interleaved(nr_cols, xparts, nr_threads, xnodes);
	vec_init_rand_range(x, (spx_value_t) -1000, (spx_value_t) 1000);

	// Allocate an interleaved y.
	y = vec_create_interleaved(nr_rows, yparts, nr_threads, ynodes);
	vec_init(y, 0);

#   if defined(SPM_NUMA) && defined(NUMA_CHECKS)
	int alloc_err = 0;
	alloc_err = check_interleaved(x->elements, xparts, nr_threads, xnodes);
	print_alloc_status("input vector", alloc_err);
	alloc_err = check_interleaved(y->elements, yparts, nr_threads, ynodes);
	print_alloc_status("output vector", alloc_err);
#   endif
#else
	x = vec_create(nr_cols, NULL);
	y = vec_create(nr_rows, NULL);
#endif

	vec_init_rand_range(x, (spx_value_t) -1000, (spx_value_t) 1000);

#ifdef DISABLE_POOL
	/*
	 * spawn two kind of threads:
	 *	- 1	: do_spmv_thread_main_swap -> computes and does swap.
	 *	- N-1 : do_spmv_thread -> just computes.
	 */
    boost::barrier cur_barrier(nr_threads);
    vector<std::shared_ptr<boost::thread> > threads(nr_threads);
    threads[0] = make_shared<boost::thread>
        (boost::bind(&do_spmv_thread_main_swap, spm_mt->spm_threads,
                     boost::ref(cur_barrier)));
    for (size_t i = 1; i < nr_threads; ++i) {
        threads[i] = make_shared<boost::thread>
            (boost::bind(&do_spmv_thread, spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
#else
    ThreadPool &pool = ThreadPool::GetInstance();
    pool.InitThreads(nr_threads - 1);
    tsc_t total_tsc;
    tsc_init(&total_tsc);
    tsc_start(&total_tsc);
    for (size_t i = 0; i < loops; i++) {
        MatVecMult(spm_mt, x, 1, y, 1);
    }
    tsc_pause(&total_tsc);
    secs = tsc_getsecs(&total_tsc);
    tsc_shut(&total_tsc);
#endif

#ifdef SPMV_PRFCNT
	// Report performance counters for every thread.
	for (size_t i = 0; i < nr_threads; i++) {
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

void spmv_check_mt(csx::CSR<uindex_t, spx_value_t> *csr, spm_mt_t *spm_mt,
                   size_t loops, size_t nr_rows, size_t nr_cols)
{
	spx_vector_t *y_tmp;
	nr_loops = loops;
    nr_threads = spm_mt->nr_threads;

#ifdef SPM_NUMA
	size_t *parts = (size_t *) malloc(sizeof(*parts)*nr_threads);
	int *nodes = (int *) malloc(sizeof(*nodes)*nr_threads);
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	for (size_t i = 0; i < spm_mt->nr_threads; i++) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
		parts[i] = spm->nr_rows * sizeof(spx_value_t);
		nodes[i] = spm->node;
	}

	int node = spm_mt->spm_threads[0].node;
	x = vec_create_onnode(nr_cols, node);
	y = vec_create_interleaved(nr_rows, parts, nr_threads, nodes);
	y_tmp = vec_create(nr_rows, NULL);
#else
	x = vec_create(nr_cols, NULL);
	y = vec_create(nr_rows, NULL);
	y_tmp = vec_create(nr_rows, NULL);
#endif
	vec_init_rand_range(x, (spx_value_t) -1000, (spx_value_t) 1000);
	vec_init(y, (spx_value_t) 0);
    vec_init(y_tmp, (spx_value_t) 0);

    vector<std::shared_ptr<boost::thread> > threads(nr_threads);
    boost::barrier cur_barrier(nr_threads + 1);

	for (size_t i = 0; i < nr_threads; ++i) {
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
        threads[i] = make_shared<boost::thread>
            (boost::bind(&do_spmv_thread, spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
	}

	for (size_t i = 0; i < loops; i++) {
		cur_barrier.wait();
		cur_barrier.wait();
		csr_spmv(csr, x, y_tmp);
		if (vec_compare(y_tmp, y) < 0)
			exit(1);
	}

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

#ifdef SPMV_PRFCNT
	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
		prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
		prfcnt_shut(prfcnt);
		free(prfcnt);
	}
#endif

	vec_destroy(x);
	vec_destroy(y);
	vec_destroy(y_tmp);
#ifdef SPM_NUMA
	free(nodes);
	free(parts);
#endif
}

float spmv_bench_sym_mt(spm_mt_t *spm_mt, size_t loops, size_t nr_rows,
                        size_t nr_cols)
{
	secs = 0.0;
	nr_loops = loops;
	nr_threads = spm_mt->nr_threads;
	n = nr_rows;

#ifdef SPMV_PRFCNT
	for (size_t i = 0; i < nr_threads; i++) {
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

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm = spm_mt->spm_threads + i;
		xparts[i] = spm->nr_rows * sizeof(spx_value_t);
		yparts[i] = spm->nr_rows * sizeof(spx_value_t);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
    }

	// Allocate an interleaved x.
	x = vec_create_interleaved(nr_rows, xparts, nr_threads, xnodes);
	vec_init_rand_range(x, (spx_value_t) -1000, (spx_value_t) 1000);

	// Allocate an interleaved y.
	y = vec_create_interleaved(nr_rows, yparts, nr_threads, ynodes);
	vec_init(y, 0);

#   if defined(SPM_NUMA) && defined(NUMA_CHECKS)
	int alloc_err = 0;
	alloc_err = check_interleaved(x->elements, xparts, nr_threads, xnodes);
	print_alloc_status("input vector", alloc_err);
	alloc_err = check_interleaved(y->elements, yparts, nr_threads, ynodes);
	print_alloc_status("output vector", alloc_err);
#   endif
#else
	x = vec_create(n, NULL);
	y = vec_create(n, NULL);
#endif

	tmp = (spx_vector_t **) malloc(nr_threads * sizeof(spx_vector_t *));
	if (!tmp) {
		perror("malloc");
		exit(1);
	}

    tmp[0] = y;
	for (size_t i = 1; i < nr_threads; i++) {
#ifdef SPM_NUMA
		int tnode = spm_mt->spm_threads[i].node;
		tmp[i] = vec_create_onnode(n, tnode);
#else
		tmp[i] = vec_create(n, NULL);
#endif
    }

	for (size_t i = 1; i < nr_threads; i++) {
		for (size_t j = 1; j < n; j++)
			tmp[i]->elements[j] = 0;
	}

#if defined(SPM_NUMA) && defined(NUMA_CHECKS)
	alloc_err = 0;
	for (size_t i = 1; i < nr_threads; i++) {
		int tnode = spm_mt->spm_threads[i].node;
		alloc_err += check_region(tmp[i]->elements, n*sizeof(spx_value_t),
                                  tnode);
	}
	print_alloc_status("local buffers", alloc_err);
#endif

#ifdef DISABLE_POOL
    boost::barrier cur_barrier(nr_threads);
    vector<std::shared_ptr<boost::thread> > threads(nr_threads);
    threads[0] = make_shared<boost::thread>
        (boost::bind(&do_spmv_thread_sym_main_swap, spm_mt->spm_threads,
                     boost::ref(cur_barrier)));
    for (size_t i = 1; i < nr_threads; ++i) {
        threads[i] = make_shared<boost::thread>
            (boost::bind(&do_spmv_thread_sym, spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
#else
    spm_mt->local_buffers =
        (spx_vector_t **) malloc(nr_threads*sizeof(spx_vector_t *));

    unsigned int i;
#   ifdef SPM_NUMA
    for (i = 1; i < nr_threads; i++) {
        int node = spm_mt->spm_threads[i].node;
        spm_mt->local_buffers[i] = vec_create_onnode(n, node);
    }
#   else
    for (i = 1; i < nr_threads; i++)
        spm_mt->local_buffers[i] = vec_create(n, NULL);
#   endif

    ThreadPool &pool = ThreadPool::GetInstance();
    pool.InitThreads(nr_threads - 1);
    tsc_t total_tsc;
    tsc_init(&total_tsc);
    tsc_start(&total_tsc);
    for (size_t i = 0; i < loops; i++) {
        MatVecMult_sym(spm_mt, x, 1, y, 1);
    }
    tsc_pause(&total_tsc);
    secs = tsc_getsecs(&total_tsc);
    tsc_shut(&total_tsc);
#endif

#ifdef SPMV_PRFCNT
	// Report performance counters for every thread.
	for (size_t i = 0; i < nr_threads; i++) {
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
	for (size_t i = 1; i < nr_threads; i++)
		vec_destroy(tmp[i]);
	free(tmp);
#ifdef SPM_NUMA
    free(xparts);
    free(xnodes);
    free(yparts);
    free(ynodes);
#endif

	return secs;
}

void spmv_check_sym_mt(csx::CSR<uindex_t, spx_value_t> *spm, spm_mt_t *spm_mt,
                       size_t loops, size_t nr_rows, size_t nr_cols)
{
	spx_vector_t *y_tmp;

	nr_loops = loops;
	nr_threads = spm_mt->nr_threads;
	n = nr_rows;
    assert(nr_rows == nr_cols);

#ifdef SPM_NUMA
	size_t *parts = (size_t *) malloc(nr_threads * sizeof(*parts));
	int *nodes = (int *) malloc(nr_threads * sizeof(*nodes));
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
		parts[i] = spm->nr_rows * sizeof(spx_value_t);
		nodes[i] = spm->node;
	}

	x = vec_create_interleaved(n, parts, nr_threads, nodes);
	y = vec_create_interleaved(n, parts, nr_threads, nodes);
	y_tmp = vec_create(n, NULL);
#else
	x = vec_create(n, NULL);
	y = vec_create(n, NULL);
	y_tmp = vec_create(n, NULL);
#endif
    vec_init_rand_range(x, (spx_value_t) -1000, (spx_value_t) 1000);
    vec_init(y, (spx_value_t) 0);
    vec_init(y_tmp, (spx_value_t) 0);

	tmp = (spx_vector_t **) malloc(nr_threads * sizeof(spx_vector_t *));
	if (!tmp) {
		perror("malloc");
		exit(1);
	}

	tmp[0] = y;
	for (size_t i = 1; i < nr_threads; i++) {
#ifdef SPM_NUMA
		int node = spm_mt->spm_threads[i].node;
        tmp[i] = vec_create_onnode(n, node);
		for (size_t j = 1; j < n; j++)
			tmp[i]->elements[j] = 0;
#else
		tmp[i] = vec_create(n, NULL);
#endif
    }

    vector<std::shared_ptr<boost::thread> > threads(nr_threads);
    boost::barrier cur_barrier(nr_threads + 1);

	for (size_t i = 0; i < nr_threads; i++) {
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
        threads[i] = make_shared<boost::thread>
            (boost::bind(&do_spmv_thread_sym, spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
	}

	for (size_t i = 0; i < loops; i++) {
		cur_barrier.wait();
		cur_barrier.wait();
		cur_barrier.wait();
		csr_spmv(spm, x, y_tmp);
		if (vec_compare(y_tmp, y) < 0)
			exit(1);
	}

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

#ifdef SPMV_PRFCNT
	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
		prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
		prfcnt_shut(prfcnt);
		free(prfcnt);
	}
#endif

	vec_destroy(x);
	vec_destroy(y);
	vec_destroy(y_tmp);
	for (size_t i = 1; i < nr_threads; i++)
		vec_destroy(tmp[i]);
	free(tmp);
#ifdef SPM_NUMA
    free(parts);
    free(nodes);
#endif
}
