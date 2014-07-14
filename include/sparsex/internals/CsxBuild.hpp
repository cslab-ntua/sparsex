/*
 * \file CsxBuild.hpp
 *
 * \brief Front-end utilities for building CSX
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_BUILD_HPP
#define SPARSEX_INTERNALS_CSX_BUILD_HPP

#include <sparsex/internals/Affinity.hpp>
#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Csx.hpp>
#include <sparsex/internals/EncodingManager.hpp>
#include <sparsex/internals/Jit.hpp>
#include <sparsex/internals/Runtime.hpp>
#include <sparsex/internals/SparseInternal.hpp>
#include <sparsex/internals/SparsePartition.hpp>
#include <sparsex/internals/SpmMt.hpp>
#include <numa.h>
#include <boost/thread/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace std;

namespace sparsex {
namespace csx {

/**
 *  Routine responsible for making a map for the symmetric representation of
 *  sparse matrices.
 *  
 *  @param spm_mt  parameters of the multithreaded execution.
 *  @param spm_sym the multithreaded CSX-Sym matrix.
 */
template<typename IndexType, typename ValueType>
void MakeMap(spm_mt_t *spm_mt,
             SparseInternal<SparsePartitionSym<IndexType, ValueType> > *spms);

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
void Compile(spm_mt_t *spm_mt, vector<ThreadContext<InternalType> > &data, 
             size_t nr_threads, bool symmetric);

/**
 *  Parallel Preprocessing.
 *
 *  @param data         thread specific data.
 *  @param csx_config   CSX configuration.
 */
template<typename IndexType, typename ValueType>
void PreprocessThreadSym(ThreadContext<SparsePartitionSym<IndexType,
                                                          ValueType> >& data);
template<typename IndexType, typename ValueType>
void PreprocessThread(ThreadContext<SparsePartition<IndexType,
                                                    ValueType> > &data);

/**
 *  Routines responsible for the construction of CSX/CSX-Sym.
 */
template<typename IndexType, typename ValueType>
void DoBuild(SparseInternal<SparsePartition<IndexType, ValueType> > *spms,
                 spm_mt *spm_mt);
template<typename IndexType, typename ValueType>
void DoBuildSym(SparseInternal<SparsePartitionSym<IndexType, ValueType> >*spms,
                spm_mt *spm_mt);
template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(SparseInternal<SparsePartition<IndexType, ValueType> > *spms);
template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsxSym(SparseInternal<SparsePartitionSym<IndexType,
                                                        ValueType> > *spms);

template<typename InternalType>
void Compile(spm_mt_t *spm_mt, vector<ThreadContext<InternalType> > &data, 
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
        data[i].GetBuffer() << "==== ENCODED PATTERNS ====\n";
        Jits[i]->GenCode(data[i].GetBuffer());    
        LOG_VERBOSE << data[i].GetBuffer().str();
    }

    for (size_t i = 0; i < nr_threads; i++)
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

    // Cleanup
    for (size_t i = 0; i < nr_threads; i++)
        delete Jits[i];
    delete[] Jits;
}

// template<typename InternalType>
// void Compile(spm_mt_t *spm_mt, vector<ThreadContext<InternalType> > &data, 
//              size_t nr_threads, bool symmetric)
// {
//     typedef typename InternalType::idx_t IndexType;
//     typedef typename InternalType::val_t ValueType;

//     CsxExecutionEngine &engine = CsxJitInit();
//     vector<CsxJit<IndexType, ValueType> > Jits;
//     for (size_t i = 0; i < nr_threads; ++i) {
//         Jits.push_back(CsxJit<IndexType, ValueType>(data[i].GetCsxManager(),
//                                                     &engine, i, symmetric));
//         Jits[i].GenCode(data[i].GetBuffer());    
//         LOG_VERBOSE << data[i].GetBuffer().str();
//     }

//     for (size_t i = 0; i < nr_threads; i++)
//         spm_mt->spm_threads[i].spmv_fn = Jits[i].GetSpmvFn();
// }

template<typename IndexType, typename ValueType>
void PreprocessThread(
    ThreadContext<SparsePartition<IndexType, ValueType> > &data)
{
    RuntimeConfiguration &rtconfig = RuntimeConfiguration::GetInstance();

    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());

    // Search for encodings
    EncodingManager<IndexType, ValueType> *DrleMg;
    data.GetBuffer() << "==> Thread: #" << data.GetId() << "\n";
    data.GetBuffer() << "==== ENCODING STATISTICS ====\n";

    // Initialize the Encoding manager.
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

    timing::Timer timer;
    timer.Start();
    CsxMatrix<IndexType, ValueType> *csx = data.GetCsxManager()->MakeCsx(false);
    timer.Pause();  // csx creation time
    double csx_create_time = timer.ElapsedTime();
    data.GetBuffer() << "CSX build: " << csx_create_time << "\n";

    data.GetSpmEncoded()->spm = csx;
    data.GetSpmEncoded()->nr_rows = csx->nrows;
    data.GetSpmEncoded()->row_start = csx->row_start;
    assert(csx->ctl);
    assert(((csx_matrix_t *)data.GetSpmEncoded()->spm)->ctl);

#if SPX_USE_NUMA
    int alloc_err = 0;
    alloc_err = check_region(csx->ctl, csx->ctl_size * sizeof(uint8_t),
                             data.GetSpmEncoded()->node);

    data.GetBuffer() << "==== ALLOCATION CHECKS ====\n";
    data.GetBuffer() << " * allocation check for ctl field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << "\n";

    alloc_err = check_region(csx->values, csx->nnz * sizeof(*csx->values),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << " * allocation check for values field... " 
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << "\n";
#endif

    delete DrleMg;
}

template<typename IndexType, typename ValueType>
void PreprocessThreadSym(
    ThreadContext<SparsePartitionSym<IndexType, ValueType> >& data)
{
    RuntimeConfiguration &rtconfig = RuntimeConfiguration::GetInstance();

    // Set CPU affinity.
    setaffinity_oncpu(data.GetCpu());

    EncodingManager<IndexType, ValueType> *DrleMg1;
    EncodingManager<IndexType, ValueType> *DrleMg2;

    data.GetBuffer() << "==> Thread: #" << data.GetId() << "\n";
    data.GetBuffer() << "==== ENCODING STATISTICS ====\n";
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
        if (data.GetId()) DrleMg1->EncodeAll(data.GetBuffer());
        DrleMg2->EncodeAll(data.GetBuffer());
    }

    data.GetSpm()->MergeMatrix();

    timing::Timer timer;
    timer.Start();
    CsxSymMatrix<IndexType, ValueType> *csx = data.GetCsxManager()->MakeCsxSym();
    timer.Pause();  // csx creation time
    double csx_create_time = timer.ElapsedTime();
    data.GetBuffer() << "CSX build: " << csx_create_time << "\n";

    data.GetSpmEncoded()->spm = csx;
    data.GetSpmEncoded()->row_start = csx->lower_matrix->row_start;
    data.GetSpmEncoded()->nr_rows = csx->lower_matrix->nrows;
   
#if SPX_USE_NUMA
    int alloc_err = check_region(csx->lower_matrix->ctl,
                                 csx->lower_matrix->ctl_size * sizeof(uint8_t),
                                 data.GetSpmEncoded()->node);
    data.GetBuffer() << "==== ALLOCATION CHECKS ====\n";
    data.GetBuffer() << " * allocation check for ctl field... "
                     << ((alloc_err) ?
                         "FAILED (see above for more info)" : "DONE")
                     << "\n";

    alloc_err = check_region(csx->lower_matrix->values,
                             csx->lower_matrix->nnz * sizeof(ValueType),
                             data.GetSpmEncoded()->node);
    data.GetBuffer() << " * allocation check for values field... "
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
void DoBuild(SparseInternal<SparsePartition<IndexType, ValueType> > *spms,
             spm_mt *spm_mt)
{
    RuntimeConfiguration &rt_config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    size_t nr_threads = rt_context.GetNrThreads();

    // Setup thread context
    vector<ThreadContext<SparsePartition
                         <IndexType, ValueType> > > mt_context(nr_threads);
    for (size_t i = 0; i < nr_threads; i++) {
        mt_context[i].SetId(i);
        mt_context[i].SetCpu(rt_context.GetAffinity(i));
        mt_context[i].SetNode(numa_node_of_cpu(rt_context.GetAffinity(i)));
        mt_context[i].SetData(spms, spm_mt); 
    }

    // Start preprocessing
    vector<boost::shared_ptr<boost::thread> > threads(nr_threads - 1);
    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i] = boost::make_shared<boost::thread>(
            PreprocessThread<IndexType, ValueType>, 
            boost::ref(mt_context[i + 1]));
    }

    // Let main thread do some work
    PreprocessThread<IndexType, ValueType>(mt_context[0]);
    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i]->join();
    }

    // CSX JIT compilation
    Compile<SparsePartition<IndexType, ValueType> >(
        spm_mt, mt_context, nr_threads,
        rt_config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric));
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsx(SparseInternal<SparsePartition<IndexType, ValueType> > *spms)
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

template<typename IndexType, typename ValueType>
void DoBuildSym(SparseInternal<SparsePartitionSym<IndexType, ValueType> > *spms_sym,
                spm_mt *spm_mt)
{
    RuntimeConfiguration &rt_config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    size_t nr_threads = rt_context.GetNrThreads();

    // Setup thread context
    vector<ThreadContext<SparsePartitionSym
                         <IndexType, ValueType> > > mt_context(nr_threads);
    for (size_t i = 0; i < nr_threads; i++) {
        mt_context[i].SetId(i);
        mt_context[i].SetCpu(rt_context.GetAffinity(i));
        mt_context[i].SetNode(numa_node_of_cpu(rt_context.GetAffinity(i)));
        mt_context[i].SetData(spms_sym, spm_mt); 
    }

    // Start preprocessing
    vector<boost::shared_ptr<boost::thread> > threads(nr_threads - 1);
    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i] = boost::make_shared<boost::thread>(
            PreprocessThreadSym<IndexType, ValueType>, 
            boost::ref(mt_context[i + 1]));
    }

    // Let main thread do some work
    PreprocessThreadSym<IndexType, ValueType>(mt_context[0]);
    for (size_t i = 0; i < nr_threads - 1; ++i) {
        threads[i]->join();
    }

    // CSX JIT compilation
    Compile<SparsePartitionSym<IndexType, ValueType> >(
        spm_mt, mt_context, nr_threads,
        rt_config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric));
}

template<typename IndexType, typename ValueType>
spm_mt_t *BuildCsxSym(SparseInternal<SparsePartitionSym<IndexType, ValueType> > *spms)
{
    spm_mt_t *spm_mt;
    RuntimeContext &rt_context = RuntimeContext::GetInstance();

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(rt_context.GetAffinity(0));
    // Initialization of the multithreaded sparse matrix representation
    spm_mt = PrepareSpmMt(); 
    // Switch Reduction Phase
    MakeMap(spm_mt, spms);
    DoBuildSym(spms, spm_mt);

    return spm_mt;
}

template<typename IndexType, typename ValueType>
void MakeMap(spm_mt_t *spm_mt,
             SparseInternal<SparsePartitionSym<IndexType, ValueType> > *spms)
{
    spm_mt_thread *spm_thread;
    SparsePartition<IndexType, ValueType> *spm;
    unsigned int *count;
    bool **initial_map;
    map_t *map;
    uint32_t start, end, col, total_count, temp_count, limit;
    uint32_t ncpus = spm_mt->nr_threads;
    uint32_t n = spms->GetPartition(0)->GetLowerMatrix()->GetNrCols();
#if SPX_USE_NUMA
    NumaAllocator &numa_alloc = NumaAllocator::GetInstance();
    int node;
#endif

    ///> Init initial_map.
    count = new unsigned int[n];
    initial_map = new bool*[ncpus];
    for (unsigned int i = 0; i < ncpus; i++)
        initial_map[i] = new bool[n];
    
    for (unsigned int i = 0; i < n; i++)
        count[i] = 0;
        
    for (unsigned int i = 0; i < ncpus; i++)
        for (unsigned int j = 0; j < n; j++)
            initial_map[i][j] = 0;

    ///> Fill initial map.
    for (unsigned int i = 0; i < ncpus; i++) {
        spm = spms->GetPartition(i)->GetLowerMatrix();

        start = spm->GetRowStart();
        end = spm->GetRowStart() + spm->GetRowptrSize() - 1;
        for (uint32_t j = 0; j < (uint32_t) spm->GetRowptrSize() - 1; j++) {
            typename SparsePartition<IndexType, ValueType>::iterator ri =
                spm->begin(j);
            typename SparsePartition<IndexType, ValueType>::iterator re =
                spm->end(j);
            for (; ri != re; ++ri) {
                col = (*ri).GetCol();
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
            cout << initial_map[i][j] << " ";
        cout << endl;
    }
    for (unsigned int i = 0; i < n; i++)
        cout << count[i] << " ";
    cout << endl;
    */
    
    ///> Make map.
    end = 0;
    for (unsigned int i = 0; i < ncpus - 1; i++) {
        spm_thread = spm_mt->spm_threads + i;
#if SPX_USE_NUMA
        node = spm_thread->node;
        map = new (numa_alloc, node) map_t;
#else
        map = new map_t;
#endif

        start = end;
        limit = total_count / (ncpus - i);
        temp_count = 0;
        while (temp_count < limit)
            temp_count += count[end++];
        total_count -= temp_count;
        map->length = temp_count;
#if SPX_USE_NUMA
        map->cpus = new (numa_alloc, node) unsigned int[temp_count];
        map->elems_pos = new (numa_alloc, node) unsigned int[temp_count];
#else
        map->cpus = new unsigned int[temp_count];
        map->elems_pos = new unsigned int[temp_count];
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

#if SPX_USE_NUMA
    node = spm_thread->node;

    map = new (numa_alloc, node) map_t;
    if (temp_count) {
        map->cpus = new (numa_alloc, node) unsigned int[temp_count];
        map->elems_pos = new (numa_alloc, node) unsigned int[temp_count];
    }
#else
    map = new map_t;
    if (temp_count) {
        map->cpus = new unsigned int[temp_count];
        map->elems_pos = new unsigned int[temp_count];
    }
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

#if SPX_USE_NUMA
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

    stringstream os;
    os << "allocation check for map... "
       << ((alloc_err) ? "FAILED (see above for more info)" : "DONE")
       << "\n";
    LOG_VERBOSE << os.str();
#endif
    ///> Print final map.  
    // for (unsigned int i = 0; i < ncpus; i++) {
    //     spm_thread = spm_mt->spm_threads + i;
    //     map = spm_thread->map;
    //     cout << "Thread " << i << endl;
    //     for (unsigned int j = 0; j < map->length; j++) {
    //         cout << "(" << map->cpus[j] << ", " << map->elems_pos[j]
    //                   << ")" << endl;
    //     }
    //     cout << endl;
    // }
    
    
    ///> Release initial map.
    for (unsigned int i = 0; i < ncpus; i++)
        delete[] initial_map[i];

    delete[] initial_map;
    delete[] count;
}    

} // end of namespace csx
} // end of namespace sparsex

#endif // SPARSEX_INTERNALS_CSX_BUILD_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
