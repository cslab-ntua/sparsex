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

#if 0   // SYM
/**
 *  Routine responsible for making a map for the symmetric representation of
 *  sparse matrices.
 *  
 *  @param spm_mt  parameters of the multithreaded execution.
 *  @param spm_sym the multithreaded CSX-Sym matrix.
 */
void MakeMap(spm_mt_t *spm_mt, SparsePartitionSym<uint64_t, double> *spm_sym);

/**
 *  Find the size of the map including the values accessed by it.
 *
 *  @param spm_mt  the complete multithreaded CSX-Sym matrix.
 *  @return        the size of the map.
 */
uint64_t MapSize(void *spm);
#endif  // SYM    

/**
 *  Parallel Preprocessing.
 *
 *  @param threadConfig parameters of matrix needed by each thread.
 *  @param csxConfig    CSX configuration
 */
template<typename IndexType, typename ValueType>
void PreprocessThread(
    ThreadContext<SparsePartition<IndexType, ValueType> > &data)
{
    RuntimeConfiguration &rtconfig = RuntimeConfiguration::GetInstance();

    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());
    
    EncodingManager<IndexType, ValueType> *DrleMg;
    data.GetBuffer() << "==> Thread: #" << data.GetId() << std::endl;

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
        cout << "hello\n";
        DrleMg->RemoveIgnore(encseq);
        DrleMg->EncodeAll(data.GetBuffer());
    }

    // DrleMg->MakeEncodeTree();
    csx_double_t *csx = data.GetCsxManager()->MakeCsx(false);
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

#if 0   // SYM
template<typename IndexType, typename ValueType>
void PreprocessThreadSym(
    ThreadContext<SparsePartitionSym<IndexType, ValueType> >& data)
{
    RuntimeConfiguration &rtconfig = RuntimeConfiguration::GetInstance();

    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());
    
    EncodingManager<IndexType, ValueType> *DrleMg1;
    EncodingManager<IndexType, ValueType> *DrleMg2;

    data.GetBuffer() << "==> Thread: #" << data.GetId() << std::endl;
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

    csx_double_sym_t *csx = data.GetCsxManager()->MakeCsxSym();
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
#endif  // SYM    

/**
 *  Routine responsible for retrieving the CSX or CSX-Sym sparse matrix format
 *  according to the command line parameters.
 *  
 *  @param mmf_fname    name of the sparse matrix file.
 *  @param spms         an initial sparse matrix format if it exists.
 *  @param spms_sym     an initial symmetric sparse matrix format if it exists.
 *  @return             the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
spm_mt_t *PrepareSpmMt();

template<typename IndexType, typename ValueType>
void DoBuild(SparseInternal<IndexType, ValueType> *spms, spm_mt *spm_mt)
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
    CsxJit<IndexType, ValueType> **Jits =
        new CsxJit<IndexType, ValueType>*[nr_threads];
    for (size_t i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit<IndexType, ValueType>(
            mt_context[i].GetCsxManager(),
            &rt_context.GetEngine(), i,
            rt_config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric));
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
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(SparseInternal<IndexType, ValueType> *spms)
{
    spm_mt_t *spm_mt;
    RuntimeContext &rt_context = RuntimeContext::GetInstance();

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_context.GetAffinity(0));
    // Initialization of the multithreaded sparse matrix representation
    spm_mt = PrepareSpmMt();
    DoBuild(spms, spm_mt);
    return spm_mt;
}

#if 0   // SYM
template<typename IndexType, typename ValueType>
void DoBuild(SparsePartitionSym<IndexType, ValueType> *spms_sym,
             spm_mt *spm_mt);
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
    CsxJit **Jits = new CsxJit*[nr_threads];
    for (size_t i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit(
            mt_context[i].GetCsxManager(),
            &rt_context.GetEngine(), i,
            rt_config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric));
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
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsxSym(SparsePartitionSym<IndexType, ValueType> *spms_sym)
{
    spm_mt_t *spm_mt;

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_context.GetAffinity(0));
    // Initialization of the multithreaded sparse matrix representation
    spm_mt = PrepareSpmMt(); 
    // Switch Reduction Phase
    MakeMap(spm_mt, spms_sym);
    DoBuild(spms_sym, spm_mt);

    return spm_mt;
}
#endif  // SYM    

/**
 *  Deallocation of CSX or CSX-Sym sparse matrix.
 *  
 *  @param spm_mt the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
void PutSpmMt(spm_mt_t *spm_mt);

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

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX form.
 *
 *  @param spm_mt  the sparse matrix in CSX format.
 */
uint64_t CsxSize(void *spm_mt);

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX-Sym form.
 *
 *  @param spm_mt  the sparse matrix in CSX-Sym format.
 */
uint64_t CsxSymSize(void *spm_mt);

#endif // SPMV_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
