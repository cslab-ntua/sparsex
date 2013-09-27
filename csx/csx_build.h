/* -*- C++ -*-
 * 
 * csx_build.h -- Front-end utilities for building CSX.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_BUILD_H__
#define CSX_BUILD_H__

#include "sparse_internal.h"
#include "sparse_util.h"
#include "sparse_partition.h"
#include "runtime.h"
#include "affinity.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"
#include "spm_mt.h"
#include "timer.h"

#include <iostream>
#include <numa.h>
#include <numaif.h>
#include <boost/thread/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

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
spm_mt_t *PrepareSpmMt();

/**
 *  Routine responsible for the Just-In-Time Compilation.
 */
template<typename InternalType>
void Compile(spm_mt_t *spm_mt, std::vector<ThreadContext<InternalType> > &data, 
             size_t nr_threads, bool symmetric);

/**
 *  Parallel Preprocessing.
 *
 *  @param data         thread specific data.
 *  @param csx_config   CSX configuration.
 */
template<typename IndexType, typename ValueType>
void PreprocessThreadSym(ThreadContext<SparsePartitionSym<IndexType, ValueType> >& data);
template<typename IndexType, typename ValueType>
void PreprocessThread(ThreadContext<SparsePartition<IndexType, ValueType> > &data);

/**
 *  Routines responsible for the construction of CSX/CSX-Sym.
 */
template<typename IndexType, typename ValueType>
double DoBuild(SparseInternal<IndexType, ValueType> *spms, spm_mt *spm_mt);
template<typename IndexType, typename ValueType>
double DoBuildSym(SparsePartitionSym<IndexType, ValueType> *spms_sym,
                  spm_mt *spm_mt);
template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(SparseInternal<IndexType, ValueType> *spms, double &time);

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsxSym(SparsePartitionSym<IndexType, ValueType> *spms_sym,
                      double &time);

template<typename InternalType>
void Compile(spm_mt_t *spm_mt, std::vector<ThreadContext<InternalType> > &data, 
             size_t nr_threads, bool symmetric)
{
    typedef typename InternalType::idx_t IndexType;
    typedef typename InternalType::val_t ValueType;

    CsxExecutionEngine &engine = CsxJitInit();
    CsxJit<IndexType, ValueType> **Jits =
        new CsxJit<IndexType, ValueType>*[nr_threads];
    for (size_t i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit<IndexType, ValueType>(data[i].GetCsxManager(),
                                                   &engine, i, symmetric);
        Jits[i]->GenCode(data[i].GetBuffer());    
        LOG_INFO << data[i].GetBuffer().str();
    }

    for (size_t i = 0; i < nr_threads; i++)
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

    // Cleanup
    for (size_t i = 0; i < nr_threads; i++)
        delete Jits[i];
    delete[] Jits;
}

template<typename IndexType, typename ValueType>
void PreprocessThread(ThreadContext<SparsePartition<IndexType, ValueType> > &data)
{
    RuntimeConfiguration &rtconfig = RuntimeConfiguration::GetInstance();

    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());
    
    EncodingManager<IndexType, ValueType> *DrleMg;
    data.GetBuffer() << "==> Thread: #" << data.GetId() << "\n";

    // Initialize the DRLE manager.
    DrleMg = new EncodingManager<IndexType, ValueType>(data.GetSpm(), rtconfig);
                                  
    // Adjust the ignore settings properly.
    EncodingSequence encseq(rtconfig.GetProperty<string>(
                                RuntimeConfiguration::PreprocXform));

    // If the user supplies the deltas choices, encode the matrix
    // with the order given in XFORM_CONF, otherwise find
    // statistical data for the types in XFORM_CONF, choose the
    // best choise, encode it and proceed likewise until there is
    // no satisfying encoding.
    if (encseq.IsExplicit()) {
        DrleMg->EncodeSerial(encseq);
    } else {
        DrleMg->RemoveIgnore(encseq);
        DrleMg->EncodeAll(data.GetBuffer());
    }

    // DrleMg->MakeEncodeTree();
    csx_t<ValueType> *csx = data.GetCsxManager()->MakeCsx(false);
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
                     << "\n";

    alloc_err = check_region(csx->values, csx->nnz * sizeof(*csx->values),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for values field... " 
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << "\n";
#endif
    delete DrleMg;
}

template<typename IndexType, typename ValueType>
void PreprocessThreadSym(ThreadContext<SparsePartitionSym<IndexType, ValueType> >& data)
{
    RuntimeConfiguration &rtconfig = RuntimeConfiguration::GetInstance();

    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());
    
    EncodingManager<IndexType, ValueType> *DrleMg1;
    EncodingManager<IndexType, ValueType> *DrleMg2;

    data.GetBuffer() << "==> Thread: #" << data.GetId() << "\n";
    data.GetSpm()->DivideMatrix();
        
    // Initialize the DRLE manager
    DrleMg1 = new EncodingManager<IndexType, ValueType>
        (data.GetSpm()->GetFirstMatrix(), rtconfig);
    DrleMg2 = new EncodingManager<IndexType, ValueType>
        (data.GetSpm()->GetSecondMatrix(), rtconfig);
                     
    EncodingSequence encseq(rtconfig.GetProperty<string>(
                                RuntimeConfiguration::PreprocXform));

    // If the user supplies the deltas choices, encode the matrix
    // with the order given in XFORM_CONF, otherwise find
    // statistical data for the types in XFORM_CONF, choose the
    // best choise, encode it and proceed likewise until there is
    // no satisfying encoding.
    if (encseq.IsExplicit()) {
        DrleMg1->EncodeSerial(encseq);
        DrleMg2->EncodeSerial(encseq);
    } else {
        DrleMg1->RemoveIgnore(encseq);
        DrleMg2->RemoveIgnore(encseq);
        DrleMg1->EncodeAll(data.GetBuffer());
        DrleMg2->EncodeAll(data.GetBuffer());
    }

    data.GetSpm()->MergeMatrix();

    csx_sym_t<ValueType> *csx = data.GetCsxManager()->MakeCsxSym();
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
                     << "\n";

    alloc_err = check_region(csx->lower_matrix->values,
                             csx->lower_matrix->nnz * sizeof(ValueType),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for values field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << "\n";
                     
    alloc_err = check_region(csx->dvalues,
	                         csx->lower_matrix->nrows * sizeof(ValueType),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << "allocation check for dvalues field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << "\n";
#endif

    delete DrleMg1;
    delete DrleMg2;
}

template<typename IndexType, typename ValueType>
double DoBuild(SparseInternal<IndexType, ValueType> *spms, spm_mt *spm_mt)
{
    RuntimeConfiguration &rt_config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    size_t nr_threads = rt_context.GetNrThreads();

    // Start timer for preprocessing
    csx::Timer timer;
    timer.Start();
    // Setup thread context
    std::vector<ThreadContext<SparsePartition
                              <IndexType, ValueType> > > mt_context(nr_threads);
    for (size_t i = 0; i < nr_threads; i++) {
        mt_context[i].SetId(i);
        mt_context[i].SetCpu(rt_context.GetAffinity(i));
        mt_context[i].SetNode(numa_node_of_cpu(rt_context.GetAffinity(i)));
        mt_context[i].SetData(spms, spm_mt); 
    }

    // Start preprocessing
    std::vector<boost::shared_ptr<boost::thread> > threads(nr_threads - 1);
    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i] = boost::make_shared<boost::thread>(
            PreprocessThread<IndexType, ValueType>, 
            boost::ref(mt_context[i + 1]));
    }
    // Let main thread do some work
    PreprocessThread<IndexType, ValueType>(mt_context[0]);
    for (size_t i = 0; i < nr_threads-1; ++i) {
        threads[i]->join();
    }

    // CSX JIT compilation
    Compile<SparsePartition<IndexType, ValueType> >(
        spm_mt, mt_context, nr_threads,
        rt_config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric));

    // Preprocessing finished; stop timer
    timer.Pause();
    pre_time = timer.ElapsedTime();
    return pre_time;
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(SparseInternal<IndexType, ValueType> *spms, double &time)
{
    spm_mt_t *spm_mt;
    RuntimeContext &rt_context = RuntimeContext::GetInstance();

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_context.GetAffinity(0));
    // Initialization of the multithreaded sparse matrix representation
    spm_mt = PrepareSpmMt();
    time = DoBuild(spms, spm_mt);
    return spm_mt;
}

template<typename IndexType, typename ValueType>
double DoBuildSym(SparsePartitionSym<IndexType, ValueType> *spms_sym, spm_mt *spm_mt)
{
    
    RuntimeConfiguration &rt_config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    size_t nr_threads = rt_context.GetNrThreads();

    // Start timer for preprocessing
    csx::Timer timer;
    timer.Start();
    // Setup thread context
    std::vector<ThreadContext<SparsePartitionSym
                              <IndexType, ValueType> > > mt_context(nr_threads);
    for (size_t i = 0; i < nr_threads; i++) {
        mt_context[i].SetId(i);
        mt_context[i].SetCpu(rt_context.GetAffinity(i));
        mt_context[i].SetNode(numa_node_of_cpu(rt_context.GetAffinity(i)));
        mt_context[i].SetDataSym(spms_sym, spm_mt); 
    }
    // Start preprocessing
    std::vector<boost::shared_ptr<boost::thread> > threads(nr_threads - 1);
    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i] = boost::make_shared<boost::thread>(
            PreprocessThreadSym<IndexType, ValueType>, 
            boost::ref(mt_context[i + 1]));
    }
    // Let main thread do some work
    PreprocessThreadSym<IndexType, ValueType>(mt_context[0]);
    for (size_t i = 0; i < nr_threads-1; ++i) {
        threads[i]->join();
    }

    // CSX JIT compilation
    Compile<SparsePartitionSym<IndexType, ValueType> >(
        spm_mt, mt_context, nr_threads,
        rt_config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric));
    // Preprocessing finished; stop timer
    timer.Pause();
    pre_time = timer.ElapsedTime();

    return pre_time;
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsxSym(SparsePartitionSym<IndexType, ValueType> *spms_sym,
                      double &time)
{
    spm_mt_t *spm_mt;
    RuntimeContext &rt_context = RuntimeContext::GetInstance();

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_context.GetAffinity(0));
    // Initialization of the multithreaded sparse matrix representation
    spm_mt = PrepareSpmMt(); 
    // Switch Reduction Phase
    MakeMap(spm_mt, spms_sym);
    time = DoBuildSym(spms_sym, spm_mt);

    return spm_mt;
}

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
    
    LOG_INFO << "allocation check for map... "
             << ((alloc_err) ? "FAILED (see above for more info)" : "DONE")
             << "\n";
>>>>>>> b66e7fb93f8284b703d2cd46a40c6347a387bd06:csx/csx_build.h

#endif
    ///> Print final map.  
    // for (unsigned int i = 0; i < ncpus; i++) {
    //     spm_thread = spm_mt->spm_threads + i;
    //     map = spm_thread->map;
    //     std::cout << "Thread " << i << std::endl;
    //     for (unsigned int j = 0; j < map->length; j++) {
    //         std::cout << "(" << map->cpus[j] << ", " << map->elems_pos[j]
    //                   << ")" << std::endl;
    //     }
    //     std::cout << std::endl;
    // }
    
    
    ///> Release initial map.
    for (unsigned int i = 0; i < ncpus; i++)
        free(initial_map[i]);
    free(initial_map);
    free(count);
}

#endif // CSX_BUILD_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
