/* -*- C++ -*-
 * 
 * spmv.h -- Front-end utilities for invoking CSX.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPMV_H__
#define SPMV_H__

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <string>
#include <numa.h>
#include <sched.h>
#include <boost/thread/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "MatrixLoading.h"
#include "runtime.h"
#include "timer.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"

#include <numa.h>
#include <numaif.h>

extern "C" {
#include <libgen.h>
#include "mt_lib.h"
#include "spmv_method.h"
#include "spm_crs.h"
#include "spm_mt.h"
#ifdef SPM_NUMA
    #include "spmv_loops_mt_numa.h"
    #include "spmv_loops_sym_mt_numa.h"
    #define SPMV_CHECK_FN spmv_double_check_mt_loop_numa
    #define SPMV_BENCH_FN spmv_double_bench_mt_loop_numa
    #define SPMV_CHECK_SYM_FN spmv_double_check_sym_mt_loop_numa
    #define SPMV_BENCH_SYM_FN spmv_double_bench_sym_mt_loop_numa
#else
    #include "spmv_loops_mt.h"
    #include "spmv_loops_sym_mt.h"
    #define SPMV_CHECK_FN spmv_double_check_mt_loop
    #define SPMV_BENCH_FN spmv_double_bench_mt_loop
    #define SPMV_CHECK_SYM_FN spmv_double_check_sym_mt_loop
    #define SPMV_BENCH_SYM_FN spmv_double_bench_sym_mt_loop
#endif // SPM_NUMA
#include "ctl_ll.h"
#include <libgen.h>
}

double pre_time;

/**
 *  Routine responsible for making a map for the symmetric representation of
 *  sparse matrices.
 *  
 *  @param spm_mt  parameters of the multithreaded execution.
 *  @param spm_sym the multithreaded CSX-Sym matrix.
 */
template<typename IndexType, typename ValueType>
void MakeMap(spm_mt_t *spm_mt,
             SparsePartitionSym<IndexType, ValueType> *spm_sym);

/**
 *  Find the size of the map including the values accessed by it.
 *
 *  @param spm_mt  the complete multithreaded CSX-Sym matrix.
 *  @return        the size of the map.
 */
uint64_t MapSize(void *spm);

/**
 *  Parallel Preprocessing.
 *
 *  @param threadConfig parameters of matrix needed by each thread.
 *  @param csxConfig    CSX configuration
 */
template<typename IndexType, typename ValueType>
void PreprocessThread(ThreadContext<SparsePartition<IndexType, ValueType> > &data,
                      const CsxContext &csx_config)
{
    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());
    
    DRLE_Manager<IndexType, ValueType> *DrleMg;
    data.GetBuffer() << "==> Thread: #" << data.GetId() << std::endl;

    // Initialize the DRLE manager.
    DrleMg = new DRLE_Manager<IndexType, ValueType>
        (data.GetSpm(), 4, 255, 0.05, csx_config.GetWindowSize(),
         DRLE_Manager<IndexType, ValueType>::SPLIT_BY_NNZ,
         csx_config.GetSamplingPortion(), csx_config.GetMaxSamples(),
         csx_config.IsSplitBlocks(), false);
                                  
    // Adjust the ignore settings properly.
    DrleMg->IgnoreAll();
    for (int i = 0; *(csx_config.GetXform()+i) != -1; ++i)
        DrleMg->RemoveIgnore(static_cast<IterOrder>
                             (*(csx_config.GetXform()+i)));

    // If the user supplies the deltas choices, encode the matrix
    // with the order given in XFORM_CONF, otherwise find
    // statistical data for the types in XFORM_CONF, choose the
    // best choise, encode it and proceed likewise until there is
    // no satisfying encoding.
    if (csx_config.GetDeltas())
        DrleMg->EncodeSerial(csx_config.GetXform(), csx_config.GetDeltas());
    else
        DrleMg->EncodeAll(data.GetBuffer());

    // DrleMg->MakeEncodeTree();
    csx_double_t *csx = (csx_double_t *) data.GetCsxManager()->MakeCsx(false);
    data.GetSpmEncoded()->spm = csx;
    data.GetSpmEncoded()->nr_rows = csx->nrows;
    data.GetSpmEncoded()->row_start = csx->row_start;

#ifdef SPM_NUMA
    int alloc_err = 0;
    alloc_err = check_region(csx->ctl, csx->ctl_size * sizeof(uint8_t),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for ctl field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << std::endl;

    alloc_err = check_region(csx->values, csx->nnz * sizeof(*csx->values),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for values field... " 
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << std::endl;
#endif
    delete DrleMg;
}

template<typename IndexType, typename ValueType>
void PreprocessThreadSym(ThreadContext<SparsePartitionSym<IndexType, ValueType> >& data,
                         const CsxContext &csx_config)
{
    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());
    
    DRLE_Manager<IndexType, ValueType> *DrleMg1;
    DRLE_Manager<IndexType, ValueType> *DrleMg2;

    data.GetBuffer() << "==> Thread: #" << data.GetId() << std::endl;

    data.GetSpm()->DivideMatrix();
        
    // Initialize the DRLE manager
    DrleMg1 = new DRLE_Manager<IndexType, ValueType>
        (data.GetSpm()->GetFirstMatrix(), 4, 255-1, 0.05,
         csx_config.GetWindowSize(),
         DRLE_Manager<IndexType, ValueType>::SPLIT_BY_NNZ,
         csx_config.GetSamplingPortion(), csx_config.GetMaxSamples(), false);

    DrleMg2 = new DRLE_Manager<IndexType, ValueType>
        (data.GetSpm()->GetSecondMatrix(), 4, 255-1, 0.05,
         csx_config.GetWindowSize(),
         DRLE_Manager<IndexType, ValueType>::SPLIT_BY_NNZ,
         csx_config.GetSamplingPortion(), csx_config.GetMaxSamples(), false);
                     
    // Adjust the ignore settings properly
    DrleMg1->IgnoreAll();
    DrleMg2->IgnoreAll();
    for (int i = 0; *(csx_config.GetXform()+i) != -1; ++i) {
        DrleMg1->
            RemoveIgnore(static_cast<IterOrder>
                         (*(csx_config.GetXform()+i)));
        DrleMg2->
            RemoveIgnore(static_cast<IterOrder>
                         (*(csx_config.GetXform()+i)));
    }
       
    // If the user supplies the deltas choices, encode the matrix with the
    // order given in XFORM_CONF, otherwise find statistical data for the 
    // types in XFORM_CONF, choose the best choise, encode it and proceed 
    // likewise until there is no satisfying encoding.
    if (csx_config.GetDeltas()) {
        DrleMg1->EncodeSerial(csx_config.GetXform(), csx_config.GetDeltas());
        DrleMg2->EncodeSerial(csx_config.GetXform(), csx_config.GetDeltas());
    } else {
        DrleMg1->EncodeAll(data.GetBuffer());
        DrleMg2->EncodeAll(data.GetBuffer());
    }
    data.GetSpm()->MergeMatrix();

    csx_double_sym_t *csx = (csx_double_sym_t *) data.GetCsxManager()->MakeCsxSym();
    data.GetSpmEncoded()->spm = csx;
    data.GetSpmEncoded()->row_start = csx->lower_matrix->row_start;
    data.GetSpmEncoded()->nr_rows = csx->lower_matrix->nrows;
   
#ifdef SPM_NUMA
    int alloc_err;
        
    alloc_err = check_region(csx->lower_matrix->ctl,
                             csx->lower_matrix->ctl_size * sizeof(uint8_t),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for ctl field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << std::endl;

    alloc_err = check_region(csx->lower_matrix->values,
                             csx->lower_matrix->nnz * sizeof(ValueType),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for values field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << std::endl;
                     
    alloc_err = check_region(csx->dvalues,
	                         csx->lower_matrix->nrows * sizeof(ValueType),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for dvalues field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << std::endl;
#endif

    delete DrleMg1;
    delete DrleMg2;
}

/**
 *  Routine responsible for retrieving the CSX or CSX-Sym sparse matrix format
 *  according to the command line parameters.
 *  
 *  @param mmf_fname    name of the sparse matrix file.
 *  @param spms         an initial sparse matrix format if it exists.
 *  @param spms_sym     an initial symmetric sparse matrix format if it exists.
 *  @param run_config   runtime configuration
 *  @param csx_config   CSX configuration
 *  @return             the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
spm_mt_t *PrepareSpmMt(const RuntimeContext &rt_config,
                       const CsxContext &csx_config);

template<typename IndexType, typename ValueType>
double DoBuild(SparseInternal<IndexType, ValueType> *spms, spm_mt *spm_mt,
               const RuntimeContext &rt_config, const CsxContext &csx_config)
{
    size_t nr_threads = rt_config.GetNrThreads();

    // Start timer for preprocessing
    csx::Timer timer;
    timer.Start();

    // Setup thread context
    std::vector<ThreadContext<SparsePartition
                              <IndexType, ValueType> > > mt_context(nr_threads);

    for (size_t i = 0; i < nr_threads; i++) {
        mt_context[i].SetId(i);
        mt_context[i].SetCpu(rt_config.GetAffinity(i));
        mt_context[i].SetNode(numa_node_of_cpu(rt_config.GetAffinity(i)));
        mt_context[i].SetData(spms, spm_mt); 
    }

    // Start preprocessing
    std::vector<boost::shared_ptr<boost::thread> > threads(nr_threads - 1);

    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i] = boost::make_shared<boost::thread>
            (PreprocessThread<IndexType, ValueType>, 
             boost::ref(mt_context[i + 1]),
             boost::ref(csx_config));
    }
    // Let main thread do some work
    PreprocessThread<IndexType, ValueType>(mt_context[0], csx_config);

    for (size_t i = 0; i < nr_threads-1; ++i) {
        threads[i]->join();
    }

    // CSX JIT compilation
    CsxJit<IndexType, ValueType> **Jits =
        new CsxJit<IndexType, ValueType>*[nr_threads];
    for (size_t i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit<IndexType, ValueType>(mt_context[i].GetCsxManager(),
                                 &rt_config.GetEngine(),
                                 i, csx_config.IsSymmetric());
        Jits[i]->GenCode(mt_context[i].GetBuffer());    
        std::cout << mt_context[i].GetBuffer().str();
    }

    for (size_t i = 0; i < nr_threads; i++)
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

    // Cleanup
    delete[] Jits;

    // Preprocessing finished; stop timer
    timer.Pause();
    pre_time = timer.ElapsedTime();

    return pre_time;
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(SparseInternal<IndexType, ValueType> *spms,
                   const RuntimeContext &rt_config, const CsxContext &csx_config,
                   double &time)
{
    spm_mt_t *spm_mt;

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_config.GetAffinity(0));
    // Initialization of the multithreaded sparse matrix representation
    spm_mt = PrepareSpmMt(rt_config, csx_config); 
    time = DoBuild(spms, spm_mt, rt_config, csx_config);

    return spm_mt;
}

template<typename IndexType, typename ValueType>
double DoBuildSym(SparsePartitionSym<IndexType, ValueType> *spms_sym,
                  spm_mt *spm_mt, const RuntimeContext &rt_config,
                  const CsxContext &csx_config)
{
    size_t nr_threads = rt_config.GetNrThreads();

    // Start timer for preprocessing
    csx::Timer timer;
    timer.Start();

    // Setup thread context
    std::vector<ThreadContext<SparsePartitionSym
                              <IndexType, ValueType> > > mt_context(nr_threads);

    for (size_t i = 0; i < nr_threads; i++) {
        mt_context[i].SetId(i);
        mt_context[i].SetCpu(rt_config.GetAffinity(i));
        mt_context[i].SetNode(numa_node_of_cpu(rt_config.GetAffinity(i)));
        mt_context[i].SetDataSym(spms_sym, spm_mt); 
    }

    // Start preprocessing
    std::vector<boost::shared_ptr<boost::thread> > threads(nr_threads - 1);

    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i] = boost::make_shared<boost::thread>
            (PreprocessThreadSym<IndexType, ValueType>, 
             boost::ref(mt_context[i + 1]),
             boost::ref(csx_config));
    }
    // Let main thread do some work
    PreprocessThreadSym<IndexType, ValueType>(mt_context[0], csx_config);

    for (size_t i = 0; i < nr_threads-1; ++i) {
        threads[i]->join();
    }

    // CSX JIT compilation
    CsxJit<IndexType, ValueType> **Jits = 
        new CsxJit<IndexType, ValueType>*[nr_threads];
    for (size_t i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit<IndexType, ValueType>(mt_context[i].GetCsxManager(),
                                                   &rt_config.GetEngine(),
                                                   i, csx_config.IsSymmetric());
        Jits[i]->GenCode(mt_context[i].GetBuffer());    
        std::cout << mt_context[i].GetBuffer().str();
    }

    for (size_t i = 0; i < nr_threads; i++)
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

    // Cleanup
    delete[] Jits;

    // Preprocessing finished; stop timer
    timer.Pause();
    pre_time = timer.ElapsedTime();

    return pre_time;
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsxSym(SparsePartitionSym<IndexType, ValueType> *spms_sym,
                      const RuntimeContext &rt_config,
                      const CsxContext &csx_config, double &time)
{
    spm_mt_t *spm_mt;

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_config.GetAffinity(0));
    // Initialization of the multithreaded sparse matrix representation
    spm_mt = PrepareSpmMt(rt_config, csx_config); 
    // Switch Reduction Phase
    MakeMap(spm_mt, spms_sym);
    time = DoBuildSym(spms_sym, spm_mt, rt_config, csx_config);

    return spm_mt;
}

/*template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(const char *mmf_fname, const RuntimeContext &rt_config,
                   const CsxContext &csx_config)
{
    spm_mt_t *spm_mt;
    SparseInternal<IndexType, ValueType> *spi = NULL;
    SparsePartitionSym<IndexType, ValueType> *spms_sym = NULL;

    // Load Matrix
    if (!csx_config.IsSymmetric()) {
        spi = LoadMMF_mt<IndexType, ValueType>(mmf_fname,
                                               rt_config.GetNrThreads());
    } else {
        spms_sym = SparsePartitionSym<IndexType, ValueType>::LoadMMF_mt(mmf_fname,
                                                  rt_config.GetNrThreads());
    }
    //spm_mt = BuildCsx(spi, spms_sym, rt_config, csx_config);

    return spm_mt;
}*/

/*template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(SparseInternal<IndexType, ValueType> *spi,
                   SparsePartitionSym<IndexType, ValueType> *spms_sym,
                   const RuntimeContext &rt_config, const CsxContext &csx_config)
{
    spm_mt_t *spm_mt;
    size_t nr_threads = rt_config.GetNrThreads();

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_config.GetAffinity(0));

    // Initialization of the multithreaded sparse matrix representation
    spm_mt = (spm_mt_t *) xmalloc(sizeof(spm_mt_t));

    spm_mt->nr_threads = nr_threads;
    spm_mt->symmetric = csx_config.IsSymmetric();
    spm_mt->spm_threads =
        (spm_mt_thread_t *) xmalloc(sizeof(spm_mt_thread_t) * nr_threads);

    for (size_t i = 0; i < nr_threads; i++) {
        spm_mt->spm_threads[i].cpu = rt_config.GetAffinity(i);
        spm_mt->spm_threads[i].node =
            numa_node_of_cpu(rt_config.GetAffinity(i));
        spm_mt->spm_threads[i].id = i;
    }

    if (spms_sym) {
        // Switch Reduction Phase
        MakeMap(spm_mt, spms_sym);
    }

    // Start timer for preprocessing
    csx::Timer timer;
    timer.Start();

    // Setup thread context
    std::vector<ThreadContext<IndexType, ValueType> > mt_context(nr_threads);

    for (size_t i = 0; i < nr_threads; i++) {
        mt_context[i].SetId(i);
        mt_context[i].SetCpu(rt_config.GetAffinity(i));
        mt_context[i].SetNode(numa_node_of_cpu(rt_config.GetAffinity(i)));
        mt_context[i].SetData(spi, spms_sym, spm_mt, csx_config.IsSymmetric()); 
    }

    // Start preprocessing
    std::vector<boost::shared_ptr<boost::thread> > threads(nr_threads - 1);

    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i] = boost::make_shared<boost::thread>
            (PreprocessThread, 
             boost::ref(mt_context[i + 1]),
             boost::ref(csx_config));
    }
    
    // Let main thread do some work
    PreprocessThread(mt_context[0], csx_config);

    for (size_t i = 0; i < nr_threads-1; ++i) {
        threads[i]->join();
    }

    // CSX JIT compilation
    CsxJit **Jits = new CsxJit*[nr_threads];
    for (size_t i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit(mt_context[i].GetCsxManager(),
                                 &rt_config.GetEngine(),
                                 i, csx_config.IsSymmetric());
        Jits[i]->GenCode(mt_context[i].GetBuffer());    
        std::cout << mt_context[i].GetBuffer().str();
    }

    for (size_t i = 0; i < nr_threads; i++)
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

    // Cleanup
    delete[] Jits;

    // Preprocessing finished; stop timer
    //timer_pause(&timer);
    //pre_time = timer_secs(&timer);
    timer.Pause();
    pre_time = timer.ElapsedTime();

    return spm_mt;
    }
*/

/**
 *  Check the CSX SpMV result against the baseline single-thread CSR
 *  implementation.
 *
 *  @param the (mulithreaded) CSX matrix.
 *  @param the MMF input file.
 */
void CheckLoop(spm_mt_t *spm_mt, char *mmf_name);

/**
 *  Run CSX SpMV and record the performance information.
 */
void BenchLoop(spm_mt_t *spm_mt, char *mmf_name);

template<typename IndexType, typename ValueType>
void MakeMap(spm_mt_t *spm_mt, SparsePartitionSym<IndexType, ValueType> *spm_sym)
{
    spm_mt_thread *spm_thread;
    SparsePartition<IndexType, ValueType> *spm;
    unsigned int *count;
    bool **initial_map;
    map_t *map;
    uint32_t start, end, col, total_count, temp_count, limit;
    uint32_t ncpus = spm_mt->nr_threads;
    uint32_t n = spm_sym->GetLowerMatrix()->GetNrCols();
#ifdef SPM_NUMA
    int node;
#endif
    
    ///> Init initial_map.
    count = (unsigned int *) xmalloc(n * sizeof(unsigned int));
    initial_map = (bool **) xmalloc(ncpus * sizeof(bool *));
    for (unsigned int i = 0; i < ncpus; i++)
        initial_map[i] = (bool *) xmalloc(n * sizeof(bool));
    
    for (unsigned int i = 0; i < n; i++)
        count[i] = 0;
        
    for (unsigned int i = 0; i < ncpus; i++)
        for (unsigned int j = 0; j < n; j++)
            initial_map[i][j] = 0;

    ///> Fill initial map.
    for (unsigned int i = 0; i < ncpus; i++) {
        spm = spm_sym[i].GetLowerMatrix();

        start = spm->GetRowStart();
        end = spm->GetRowStart() + spm->GetRowptrSize() - 1;
        
        for (uint32_t j = 0; j < (uint32_t) spm->GetRowptrSize() - 1; j++) {
            for (Elem<IndexType, ValueType> *elem = spm->RowBegin(j); elem != spm->RowEnd(j);
                 elem++) {
                col = elem->col;
                assert(col < end);
                if (col < start + 1 && !initial_map[i][col]) {
                    initial_map[i][col] = 1;
                    count[col]++;
                }
            }
        }
    }
    total_count = 0;
    for (unsigned int i = 0; i < n; i++)
        total_count += count[i];
    
    ///> Print initial map.
    /*
    for (unsigned int i = 0; i < ncpus; i++) {
        for (unsigned int j = 0; j < n; j++)
            std::cout << initial_map[i][j] << " ";
        std::cout << std::endl;
    }
    for (unsigned int i = 0; i < n; i++)
        std::cout << count[i] << " ";
    std::cout << std::endl;
    */
    
    ///> Make map.
    end = 0;
    for (unsigned int i = 0; i < ncpus - 1; i++) {
        spm_thread = spm_mt->spm_threads + i;
#ifdef SPM_NUMA
        node = spm_thread->node;
        map = (map_t *) numa_alloc_onnode(sizeof(map_t), node);
#else
        map = (map_t *) xmalloc(sizeof(map_t));
#endif

        start = end;
        limit = total_count / (ncpus - i);
        temp_count = 0;
        while (temp_count < limit)
            temp_count += count[end++];
        total_count -= temp_count;
        map->length = temp_count;
#ifdef SPM_NUMA
        map->cpus = (unsigned int *) 
                        numa_alloc_onnode(temp_count * sizeof(unsigned int),
                                          node);
        map->elems_pos = (unsigned int *)
                             numa_alloc_onnode(temp_count *
                                               sizeof(unsigned int), node);
#else
        map->cpus = (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
        map->elems_pos = 
            (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
#endif        
        temp_count = 0;
        for (unsigned int j = start; j < end; j++) {
            for (unsigned int k = 0; k < ncpus; k++) {
                if (initial_map[k][j] == 1) {
                     map->cpus[temp_count] = k;
                     map->elems_pos[temp_count++] = j - 1;
                }
            }
        }
        assert(temp_count == map->length);
        
        spm_thread->map = map;
    }
    start = end;
    end = n;
    temp_count = total_count;
    
    spm_thread = spm_mt->spm_threads + ncpus - 1;

#ifdef SPM_NUMA
    node = spm_thread->node;

    map = (map_t *) numa_alloc_onnode(sizeof(map_t), node);
    map->cpus = (unsigned int *)
                    numa_alloc_onnode(temp_count * sizeof(unsigned int), node);
    map->elems_pos = (unsigned int *)
                         numa_alloc_onnode(temp_count * sizeof(unsigned int),
                                           node);
#else
    map = (map_t *) xmalloc(sizeof(map_t));
    map->cpus = (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
    map->elems_pos = (unsigned int *)
                         xmalloc(temp_count * sizeof(unsigned int));
#endif

    map->length = temp_count;
    temp_count = 0;
    for (unsigned int j = start; j < end; j++) {
        for (unsigned int k = 0; k < ncpus; k++) {
            if (initial_map[k][j] == 1) {
                map->cpus[temp_count] = k;
                map->elems_pos[temp_count++] = j - 1;
            }
        }
    }
    assert(temp_count == map->length); 
    spm_thread->map = map;
#ifdef SPM_NUMA
    int alloc_err = 0;
    for (unsigned int i = 0; i < ncpus; i++) {
        node = spm_mt->spm_threads[i].node;
        map = spm_mt->spm_threads[i].map;
        
        alloc_err += check_region(map, sizeof(map_t), node);
        alloc_err += check_region(map->cpus, 
                                    map->length * sizeof(unsigned int), node);
        alloc_err += check_region(map->elems_pos,
                                    map->length * sizeof(unsigned int), node);
    }
    
    std::cout << "allocation check for map... "
              << ((alloc_err) ?
                 "FAILED (see above for more info)" : "DONE")
              << std::endl;

#endif
    ///> Print final map.
    /*
    for (unsigned int i = 0; i < ncpus; i++) {
        spm_thread = spm_mt->spm_threads + i;
        map = spm_thread->map;
        std::cout << "Thread " << i << std::endl;
        for (unsigned int j = 0; j < map->length; j++) {
            std::cout << "(" << map->cpus[j] << ", " << map->elems_pos[j]
                      << ")" << std::endl;
        }
        std::cout << std::endl;
    }
    */
    
    ///> Release initial map.
    for (unsigned int i = 0; i < ncpus; i++)
        free(initial_map[i]);
    free(initial_map);
    free(count);
}

#endif // SPMV_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
