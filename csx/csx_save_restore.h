/* -*- C++ -*-
 *
 * csx_save_restore.h --  Saving/Restoring Csx to/from archive  
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSXSAVERESTORE_H__
#define CSXSAVERESTORE_H__

#include "csx.h"
#include "jit.h"
#include "logger.hpp"

#ifdef SPM_NUMA
#   include "affinity.h"
#   include "numa_util.h"
#   include <numa.h>
#endif

#include <boost/archive/tmpdir.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/binary_object.hpp>
#include <fstream>

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
void SaveCsx(void *spm, const char *filename, IndexType *permutation)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    csx_t<ValueType> *csx = 0;
    csx_sym_t<ValueType> *csx_sym = 0;
    ofstream file;
    bool reordered = false;

    if (!filename) {
        string default_filename(boost::archive::tmpdir());
        default_filename += "/csx_file";
        filename = default_filename.c_str();
    }
    try {
        file.open(filename, ios::binary | ios::out);
        if (!file.good()) {
            throw ios_base::failure("");
        }
        boost::archive::binary_oarchive oa(file);
        oa << spm_mt->nr_threads;
        oa << spm_mt->symmetric;
        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            if (spm_mt->symmetric) {
                csx_sym = (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;
                csx = (csx_t<ValueType> *) csx_sym->lower_matrix;
            } else {
                csx = (csx_t<ValueType> *) spm_mt->spm_threads[i].spm;
            }
            oa << spm_mt->spm_threads[i].cpu & spm_mt->spm_threads[i].id 
                & spm_mt->spm_threads[i].node;
            oa << csx->nnz & csx->ncols & csx->nrows & csx->ctl_size
                & csx->row_start;
            oa << boost::serialization::make_array(csx->values, csx->nnz);
            oa << boost::serialization::make_array(csx->ctl, csx->ctl_size);
            oa << csx->id_map;
            oa << csx->row_jumps;
            oa << boost::serialization::make_array(csx->rows_info, csx->nrows);
            if (spm_mt->symmetric) {
                oa << boost::serialization::make_array(csx_sym->dvalues, csx->nrows);
                map_t *map = spm_mt->spm_threads[i].map;
                oa << map->length;
                oa << boost::serialization::make_array(map->cpus, map->length);
                oa << boost::serialization::make_array(map->elems_pos, map->length);
            }
        }
        if (permutation != 0) {
            reordered = true;
            oa << reordered;
            oa << boost::serialization::make_array(permutation, csx->ncols);
        } else {
            oa << reordered;
        }
    } catch (ios_base::failure &e) {
        LOG_ERROR << "CSX file error\n";
        exit(1);
    } catch (boost::archive::archive_exception &e) {
        LOG_ERROR << "dumping CSX to file failed: " << e.what() << "\n";
        exit(1);
    }
    //file.close(); //Any open file is automatically closed when the ofstream is destroyed
}

template<typename IndexType, typename ValueType>
spm_mt_t *RestoreCsx(const char *filename, IndexType **permutation)
{
    unsigned int nr_threads;
    bool symmetric, reordered, full_column_indices = false;
    ifstream file;
    spm_mt_t *spm_mt = 0;
    csx_sym_t<ValueType> *csx_sym = 0;
    csx_t<ValueType> *csx = 0;

    if (!filename) {
        string default_filename(boost::archive::tmpdir());
        default_filename += "/csx_file";
        filename = default_filename.c_str();
    }
    try {
        file.open(filename, ios::binary | ios::in);
        if (!file.good()) {
            throw ios_base::failure("");
        }
        boost::archive::binary_iarchive ia(file);

        ia >> nr_threads;
        ia >> symmetric;
        // Construct spm_mt
        spm_mt = (spm_mt_t *) xmalloc(sizeof(spm_mt_t));
        spm_mt->nr_threads = nr_threads;
        spm_mt->symmetric = symmetric;
        spm_mt->spm_threads = (spm_mt_thread_t *) xmalloc
            (sizeof(spm_mt_thread_t) * nr_threads);

#ifdef SPM_NUMA
        full_column_indices = true;
        int current_node = 0;
#endif

        for (unsigned int i = 0; i < nr_threads; i++) {
            ia >> spm_mt->spm_threads[i].cpu & spm_mt->spm_threads[i].id 
                & spm_mt->spm_threads[i].node;
#ifdef SPM_NUMA
            int node = spm_mt->spm_threads[i].node;
            if (node != current_node) {
                current_node = node;
                setaffinity_oncpu(spm_mt->spm_threads[i].cpu);
            }
            if (symmetric) {
                csx_sym = (csx_sym_t<ValueType> *) numa_alloc_onnode
                    (sizeof(csx_sym_t<ValueType>), node);
            }
            csx = (csx_t<ValueType> *) alloc_onnode
                (sizeof(csx_t<ValueType>), node);
#else
            if (symmetric) {
                csx_sym = (csx_sym_t<ValueType> *) malloc
                    (sizeof(csx_sym_t<ValueType>));
            }
            csx = (csx_t<ValueType> *) malloc
                (sizeof(csx_t<ValueType>));
#endif
            if (!csx) {
                LOG_ERROR << "CSX malloc failed\n";
                exit(1);
            }
            ia >> csx->nnz & csx->ncols & csx->nrows & csx->ctl_size
                & csx->row_start;
#ifdef SPM_NUMA
            if (symmetric) {
                csx_sym->dvalues = (ValueType *) numa_alloc_onnode
                    (csx->nrows*sizeof(ValueType), node);
            }
            csx->values = (ValueType *) alloc_onnode(sizeof(ValueType)
                                                     *csx->nnz, node);
            int alloc_err = 0;
            alloc_err = check_region(csx->values, csx->nnz*sizeof(*csx->values),
                                     node);
            print_alloc_status("values", alloc_err);
            csx->ctl = (uint8_t *) alloc_onnode(sizeof(uint8_t)*csx->ctl_size,
                                                node);
            alloc_err = check_region(csx->ctl, csx->ctl_size*sizeof(uint8_t),
                                     node);
            print_alloc_status("ctl", alloc_err);
            csx->rows_info = (row_info_t *) 
                alloc_onnode(sizeof(row_info_t)*csx->nrows, node);
#else
            if (symmetric) {
                csx_sym->dvalues = (ValueType *) malloc(csx->nrows*
                                                        sizeof(ValueType));
            }
            csx->values = (ValueType *) malloc(sizeof(ValueType)*csx->nnz);
            csx->ctl = (uint8_t *) malloc(sizeof(uint8_t)*csx->ctl_size);
            csx->rows_info = (row_info_t *) 
                malloc(sizeof(row_info_t)*csx->nrows);
#endif    
            if (!csx->values || !csx->ctl) {
                LOG_ERROR << "CSX malloc failed\n";
                exit(1);
            }
            ia >> boost::serialization::make_array(csx->values, csx->nnz);
            ia >> boost::serialization::make_array(csx->ctl, csx->ctl_size);
            ia >> csx->id_map;
            ia >> csx->row_jumps;
            ia >> boost::serialization::make_array(csx->rows_info, csx->nrows);
            if (symmetric) {
                ia >> boost::serialization::make_array(csx_sym->dvalues, csx->nrows);
                map_t *map = 0;
                unsigned int length;
                ia >> length;
#ifdef SPM_NUMA
                map = (map_t *) numa_alloc_onnode(sizeof(map_t), node);
                map->cpus = (unsigned int *) 
                    numa_alloc_onnode(length * sizeof(unsigned int), node);
                map->elems_pos = (unsigned int *)
                    numa_alloc_onnode(length * sizeof(unsigned int), node);
#else
                map = (map_t *) xmalloc(sizeof(map_t));
                map->cpus = (unsigned int *) xmalloc(length * sizeof(unsigned int));
                map->elems_pos = 
                    (unsigned int *) xmalloc(length * sizeof(unsigned int));
#endif
                map->length = length;
                ia >> boost::serialization::make_array(map->cpus, length);
                ia >> boost::serialization::make_array(map->elems_pos, length);
                spm_mt->spm_threads[i].map = map;
                csx_sym->lower_matrix = csx;
                spm_mt->spm_threads[i].spm = csx_sym;
            } else {
                spm_mt->spm_threads[i].spm = csx;
            }
            spm_mt->spm_threads[i].row_start = csx->row_start;
            spm_mt->spm_threads[i].nr_rows = csx->nrows;
        }
        ia >> reordered;
        if (reordered) {
            *permutation = (IndexType *) malloc(sizeof(IndexType)*csx->ncols);
            ia >> boost::serialization::make_array(*permutation, csx->ncols);
        }

        // Initialize the CSX JIT execution engine
        CsxExecutionEngine &engine = CsxJitInit();
        CsxJit<IndexType, ValueType> **Jits =
            new CsxJit<IndexType, ValueType>*[nr_threads];
        for (size_t i = 0; i < nr_threads; ++i) {
            if (symmetric) {
                csx_sym = (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;
                csx = (csx_t<ValueType> *) csx_sym->lower_matrix;
            } else { 
                csx = (csx_t<ValueType> *) spm_mt->spm_threads[i].spm;
            }
            bool row_jumps = csx->row_jumps != 0;
            Jits[i] = new CsxJit<IndexType, ValueType>(csx, &engine,
                                                       i, symmetric, row_jumps,
                                                       full_column_indices);
            Jits[i]->GenCode(std::cout);
        }

        for (size_t i = 0; i < nr_threads; i++)
            spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

        // Cleanup
        delete[] Jits;
    } catch (ios_base::failure &e) {
        LOG_ERROR << "CSX file error\n";
        exit(1);
    } catch (boost::archive::archive_exception &e) {
        LOG_ERROR << "loading CSX from file failed: " << e.what() << "\n";
        exit(1);
    } catch (bad_alloc &e) {
        LOG_ERROR << "loading CSX from file failed: " << e.what() << "\n";
        exit(1);
    }
    //file.close(); //Any open file is automatically closed when the ifstream is destroyed

    return spm_mt;
}

#endif  // CSXSAVERESTORE_H__
