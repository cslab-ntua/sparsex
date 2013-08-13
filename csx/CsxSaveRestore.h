/* -*- C++ -*-
 *
 * CsxSaveRestore.h --  Saving/Restoring Csx to/from archive  
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSXSAVERESTORE_H__
#define CSXSAVERESTORE_H__

#include <boost/archive/tmpdir.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/binary_object.hpp>
#include <fstream>

#include "csx.h"
#include "jit.h"

#ifdef SPM_NUMA
#   include <numa.h>
#   include "numa_util.h"
#endif

using namespace std;

// malloc wrapper
#define xmalloc(x)                         \
({                                         \
    void *ret_;                            \
    ret_ = malloc(x);                      \
    if (ret_ == NULL){                     \
        std::cerr << __FUNCTION__          \
                  << " " << __FILE__       \
                  << ":" << __LINE__       \
                  << ": malloc failed\n";  \
        exit(1);                           \
    }                                      \
    ret_;                                  \
})

// Non-intrusive version
namespace boost { 
namespace serialization {

template<class Archive>
void save(Archive &ar, const row_info_t &row_info, const unsigned int version)
{
    ar & row_info.rowptr & row_info.valptr & row_info.span;
}

template<class Archive>
void load(Archive &ar, row_info_t &row_info, const unsigned int version)
{
    ar & row_info.rowptr & row_info.valptr & row_info.span;
}

template<class Archive>
inline void serialize(Archive &ar, row_info_t &row_info,
                      const unsigned int version)
{
    boost::serialization::split_free(ar, row_info, version);
}

}
}

template<typename IndexType, typename ValueType>
void SaveCsx(void *spm, const char *filename, bool symmetric,
             IndexType *permutation)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    csx_t<IndexType, ValueType> *csx = 0;
    ofstream file;
    bool reordered = false;

    if (!filename) {
        string default_filename(boost::archive::tmpdir());
        default_filename += "/csx_file";
        file.open(default_filename.c_str(), ios::binary | ios::out);
    } else {
        file.open(filename, ios::binary | ios::out);
    }

    if (file.good()) {
        boost::archive::binary_oarchive oa(file);
        oa << spm_mt->nr_threads;
        oa << symmetric;
        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            csx = (csx_t<IndexType, ValueType> *) spm_mt->spm_threads[i].spm;
            oa << csx->nnz & csx->ncols & csx->nrows & csx->ctl_size
                & csx->row_start;
            oa << boost::serialization::make_array(csx->values, csx->nnz);
            oa << boost::serialization::make_array(csx->ctl, csx->ctl_size);
            oa << csx->id_map;
            oa << csx->row_jumps;
            oa << boost::serialization::make_array(csx->rows_info, csx->nrows);
        }
        if (permutation != 0) {
            reordered = true;
            oa << reordered;
            oa << boost::serialization::make_array(permutation, csx->ncols);
        } else {
            oa << reordered;
        }
        csx = 0;
        file.close();
    } else {
        file.close();
        cerr << "Csx file error!\n";
        exit(1);
    }
}

template<typename IndexType, typename ValueType>
spm_mt_t *RestoreCsx(const char *filename, const RuntimeContext &rt_config,
                     IndexType **permutation)
{
    unsigned int nr_threads;
    bool symmetric, reordered, full_column_indices = false;
    ifstream file;
    spm_mt_t *spm_mt = 0;
    csx_t<IndexType, ValueType> *csx = 0;

    if (!filename) {
        string default_filename(boost::archive::tmpdir());
        default_filename += "/csx_file";
        file.open(default_filename.c_str(), ios::binary | ios::in);
    } else {
        file.open(filename, ios::binary | ios::in);
    }

    if (file.good()) {
        boost::archive::binary_iarchive ia(file);

        ia >> nr_threads;
        ia >> symmetric;
        // Construct spm_mt
        spm_mt = (spm_mt_t *) xmalloc(sizeof(spm_mt_t));
        spm_mt->nr_threads = nr_threads;
        spm_mt->symmetric = symmetric;
        spm_mt->spm_threads = (spm_mt_thread_t *) xmalloc
            (sizeof(spm_mt_thread_t) * nr_threads);

        for (unsigned int i = 0; i < nr_threads; i++) {
            spm_mt->spm_threads[i].cpu = rt_config.GetAffinity(i);
            spm_mt->spm_threads[i].node =
                numa_node_of_cpu(rt_config.GetAffinity(i));
            spm_mt->spm_threads[i].id = i;
        }

#ifdef SPM_NUMA
        full_column_indices = true;
#endif

        for (unsigned int i = 0; i < nr_threads; i++) {
#ifdef SPM_NUMA
            int node = numa_node_of_cpu(rt_config.GetAffinity(i));
            if (node < 0) {
                perror("numa_node_of_cpu() failed");
                exit(1);
            }
            csx = (csx_t<IndexType, ValueType> *) alloc_onnode
                (sizeof(csx_t<IndexType, ValueType>), node);
#else
            csx = (csx_t<IndexType, ValueType> *) malloc
                (sizeof(csx_t<IndexType, ValueType>));
#endif
            if (!csx) {
                std::cerr << "csx malloc failed\n";
                exit(1);
            }
            ia >> csx->nnz & csx->ncols & csx->nrows & csx->ctl_size
                & csx->row_start;
#ifdef SPM_NUMA
            csx->values = (ValueType *) alloc_onnode(sizeof(ValueType)
                                                     *csx->nnz, node);
            csx->ctl = (uint8_t *) alloc_onnode(sizeof(uint8_t)*csx->ctl_size,
                                                node);
            csx->rows_info = (row_info_t *) 
                alloc_onnode(sizeof(row_info_t)*csx->nrows, node);
#else
            csx->values = (ValueType *) malloc(sizeof(ValueType)*csx->nnz);
            csx->ctl = (uint8_t *) malloc(sizeof(uint8_t)*csx->ctl_size);
            csx->rows_info = (row_info_t *) 
                malloc(sizeof(row_info_t)*csx->nrows);
#endif    
            if (!csx->values || !csx->ctl) {
                std::cerr << "csx malloc failed\n";
                exit(1);
            }
            ia >> boost::serialization::make_array(csx->values, csx->nnz);
            ia >> boost::serialization::make_array(csx->ctl, csx->ctl_size);
            ia >> csx->id_map;
            ia >> csx->row_jumps;
            ia >> boost::serialization::make_array(csx->rows_info, csx->nrows);
            spm_mt->spm_threads[i].spm = csx;
            spm_mt->spm_threads[i].row_start = csx->row_start;
            spm_mt->spm_threads[i].nr_rows = csx->nrows;
        }
        ia >> reordered;
        if (reordered) {
            *permutation = (IndexType *) malloc(sizeof(IndexType)*csx->ncols);
            ia >> boost::serialization::make_array(*permutation, csx->ncols);
        }
        csx = 0;
        
        // CSX JIT compilation
        CsxJit<IndexType, ValueType> **Jits =
            new CsxJit<IndexType, ValueType>*[nr_threads];
        for (size_t i = 0; i < nr_threads; ++i) {
            csx_t<IndexType, ValueType> *csx = (csx_t<IndexType, ValueType> *)
                spm_mt->spm_threads[i].spm;
            bool row_jumps = csx->row_jumps != 0;
            Jits[i] = new CsxJit<IndexType, ValueType>(csx, &rt_config.GetEngine(),
                                                       i, symmetric, row_jumps,
                                                       full_column_indices);
            Jits[i]->GenCode(std::cout);
        }

        for (size_t i = 0; i < nr_threads; i++)
            spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

        // Cleanup
        delete[] Jits;

        return spm_mt;
    } else {
        cerr << "Csx file error!\n";
        exit(1);
    }
}

#endif  // CSXSAVERESTORE_H__
