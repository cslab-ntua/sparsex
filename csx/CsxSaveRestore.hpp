/*
 * CsxSaveRestore.hpp -- Saving/Restoring Csx to/from archive.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SAVE_RESTORE_HPP
#define CSX_SAVE_RESTORE_HPP

#include "Affinity.hpp"
#include "Allocators.hpp"
#include "Csx.hpp"
#include "Jit.hpp"
#include "Logger.hpp"

#ifdef SPM_NUMA
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

}   // end of namespace serialization
}   // end of namespace boost

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

        oa << spm_mt->nr_threads & spm_mt->symmetric;

        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            oa << spm_mt->spm_threads[i].cpu & spm_mt->spm_threads[i].id 
                & spm_mt->spm_threads[i].node;
#ifdef SPM_NUMA
            if (spm_mt->symmetric) {
                csx_sym = (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;
                csx = (csx_t<ValueType> *) csx_sym->lower_matrix;
            } else {
                csx = (csx_t<ValueType> *) spm_mt->spm_threads[i].spm;
            }
            oa << csx->nnz & csx->ctl_size;
#endif
        }

        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            if (spm_mt->symmetric) {
                csx_sym = (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;
                csx = (csx_t<ValueType> *) csx_sym->lower_matrix;
            } else {
                csx = (csx_t<ValueType> *) spm_mt->spm_threads[i].spm;
            }

            oa << csx->nnz & csx->ncols & csx->nrows & csx->ctl_size
                & csx->row_start;
            oa << boost::serialization::make_array(csx->values, csx->nnz);
            oa << boost::serialization::make_array(csx->ctl, csx->ctl_size);
            oa << csx->id_map & csx->row_jumps;
            oa << boost::serialization::make_array(csx->rows_info, csx->nrows);

            if (spm_mt->symmetric) {
                oa << boost::serialization::make_array(csx_sym->dvalues,
                                                       csx->nrows);
                map_t *map = spm_mt->spm_threads[i].map;
                oa << map->length;
                oa << boost::serialization::make_array(map->cpus, map->length);
                oa << boost::serialization::make_array(map->elems_pos,
                                                       map->length);
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
        LOG_ERROR << "CSX file error: " << e.what() << "\n";
        exit(1);
    } catch (boost::archive::archive_exception &e) {
        LOG_ERROR << "dumping CSX to file failed: " << e.what() << "\n";
        exit(1);
    }
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

#ifdef SPM_NUMA
    NumaAllocator &numa_alloc = NumaAllocator::GetInstance();
#endif

    setaffinity_oncpu(0);
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
        ia >> nr_threads & symmetric;

        // Construct spm_mt
        spm_mt = new spm_mt_t;
        spm_mt->nr_threads = nr_threads;
        spm_mt->symmetric = symmetric;
        spm_mt->spm_threads = new spm_mt_thread_t[nr_threads];

#ifdef SPM_NUMA
        full_column_indices = true;

        // Allocate structures with alloc_interleaved
        // std::vector<size_t> ctl_parts(nr_threads);
        // std::vector<size_t> val_parts(nr_threads);
        // std::vector<int> ctl_nodes(nr_threads);
        // std::vector<int> val_nodes(nr_threads);
        size_t *ctl_parts = new size_t[nr_threads];
        size_t *val_parts = new size_t[nr_threads];
        int *ctl_nodes = new int[nr_threads];
        int *val_nodes = new int[nr_threads];

        size_t total_nnz = 0, total_ctlsize = 0;
        size_t nnz, ctl_size;
#endif
        for (unsigned int i = 0; i < nr_threads; i++) {
            ia >> spm_mt->spm_threads[i].cpu & spm_mt->spm_threads[i].id 
                & spm_mt->spm_threads[i].node;
#ifdef SPM_NUMA
            ia >> nnz & ctl_size;
            ctl_parts[i] = ctl_size * sizeof(uint8_t);
            val_parts[i] = nnz * sizeof(ValueType);
            ctl_nodes[i] = val_nodes[i] = spm_mt->spm_threads[i].node;
            total_ctlsize += ctl_size;
            total_nnz += nnz;
#endif
        }

#ifdef SPM_NUMA
        spm_mt->interleaved = true;
        // uint8_t *ctl_interleaved = 
        //     new (numa_alloc, ctl_parts, ctl_nodes) uint8_t[total_ctlsize];      
        // ValueType *values_interleaved = 
        //     new (numa_alloc, val_parts, val_nodes) ValueType[total_nnz];
        uint8_t *ctl_interleaved = (uint8_t *) alloc_interleaved
            (total_ctlsize * sizeof(uint8_t), ctl_parts, nr_threads, ctl_nodes);
        ValueType *values_interleaved = (ValueType *) alloc_interleaved
            (total_nnz * sizeof(ValueType), val_parts, nr_threads, val_nodes);
        size_t values_index = 0, ctl_index = 0;
		int cpu = sched_getcpu();
        int node = numa_node_of_cpu(cpu);
#endif

        for (unsigned int i = 0; i < nr_threads; i++) {
#ifdef SPM_NUMA
            if (symmetric)
                csx_sym = new (numa_alloc, node) csx_sym_t<ValueType>;
            csx = new (numa_alloc, node) csx_t<ValueType>;
#else
            if (symmetric)
                csx_sym = new csx_sym_t<ValueType>;
            csx = new csx_t<ValueType>;
#endif
            ia >> csx->nnz & csx->ncols & csx->nrows & csx->ctl_size &
                csx->row_start;
#ifdef SPM_NUMA
            if (symmetric)
                csx_sym->dvalues = new (numa_alloc, node) ValueType[csx->nrows];
            csx->values = values_interleaved + values_index;
            values_index += csx->nnz;
            csx->ctl = ctl_interleaved + ctl_index;
            ctl_index += csx->ctl_size;
            csx->rows_info = new (numa_alloc, node) row_info_t[csx->nrows];
#else
            if (symmetric)
                csx_sym->dvalues = new ValueType[csx->nrows];
            csx->values = new ValueType[csx->nnz];
            csx->ctl = new uint8_t[csx->ctl_size];
            csx->rows_info = new row_info_t[csx->nrows];
#endif    
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
                map = new (numa_alloc, node) map_t;
                map->cpus = new (numa_alloc, node) unsigned int[length];
                map->elems_pos = new (numa_alloc, node) unsigned int[length];
#else
                map = new map_t;
                map->cpus = new unsigned int[length];
                map->elems_pos = new unsigned int[length];
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
            *permutation = new IndexType[csx->ncols];
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
                                                       i, symmetric, 
                                                       full_column_indices,
                                                       row_jumps);
            Jits[i]->GenCode(std::cout);
        }

        for (size_t i = 0; i < nr_threads; i++)
            spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

        // Cleanup
        for (size_t i = 0; i < nr_threads; i++)
            delete Jits[i];
        delete[] Jits;

#ifdef SPM_NUMA
#ifdef NUMA_CHECKS
        int alloc_err = check_interleaved(ctl_interleaved, ctl_parts, nr_threads,
                                          ctl_nodes);
        print_alloc_status("ctl_interleaved", alloc_err);
        alloc_err = check_interleaved(values_interleaved, val_parts, nr_threads,
				      val_nodes);
        print_alloc_status("values_interleaved", alloc_err);
#endif
        free(val_parts);
        free(ctl_parts);
        free(val_nodes);
        free(ctl_nodes);
#endif
    } catch (ios_base::failure &e) {
        LOG_ERROR << "CSX file error: " << e.what() << "\n";
        exit(1);
    } catch (boost::archive::archive_exception &e) {
        LOG_ERROR << "loading CSX from file failed: " << e.what() << "\n";
        exit(1);
    } catch (bad_alloc &e) {
        LOG_ERROR << "loading CSX from file failed: " << e.what() << "\n";
        exit(1);
    }

    return spm_mt;
}

#endif  // CSX_SAVE_RESTORE_HPP

// template<typename IndexType, typename ValueType>
// void SaveCsx(void *spm, const char *filename, IndexType *permutation)
// {
//     spm_mt_t *spm_mt = (spm_mt_t *) spm;
//     csx_t<ValueType> *csx = 0;
//     csx_sym_t<ValueType> *csx_sym = 0;
//     ofstream file;
//     bool reordered = false;

//     if (!filename) {
//         string default_filename(boost::archive::tmpdir());
//         default_filename += "/csx_file";
//         filename = default_filename.c_str();
//     }

//     try {
//         file.open(filename, ios::binary | ios::out);
//         if (!file.good()) {
//             throw ios_base::failure("");
//         }
//         boost::archive::binary_oarchive oa(file);
        
//         oa << spm_mt->nr_threads & spm_mt->symmetric;

//         for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
//             if (spm_mt->symmetric) {
//                 csx_sym = (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;
//                 csx = (csx_t<ValueType> *) csx_sym->lower_matrix;
//             } else {
//                 csx = (csx_t<ValueType> *) spm_mt->spm_threads[i].spm;
//             }
//             oa << spm_mt->spm_threads[i].cpu & spm_mt->spm_threads[i].id 
//                 & spm_mt->spm_threads[i].node;
//             oa << csx->nnz & csx->ncols & csx->nrows & csx->ctl_size
//                 & csx->row_start;
//             oa << boost::serialization::make_array(csx->values, csx->nnz);
//             oa << boost::serialization::make_array(csx->ctl, csx->ctl_size);
//             oa << csx->id_map;
//             oa << csx->row_jumps;
//             oa << boost::serialization::make_array(csx->rows_info, csx->nrows);
//             if (spm_mt->symmetric) {
//                 oa << boost::serialization::make_array(csx_sym->dvalues, csx->nrows);
//                 map_t *map = spm_mt->spm_threads[i].map;
//                 oa << map->length;
//                 oa << boost::serialization::make_array(map->cpus, map->length);
//                 oa << boost::serialization::make_array(map->elems_pos, map->length);
//             }
//         }
//         if (permutation != 0) {
//             reordered = true;
//             oa << reordered;
//             oa << boost::serialization::make_array(permutation, csx->ncols);
//         } else {
//             oa << reordered;
//         }
//     } catch (ios_base::failure &e) {
//         LOG_ERROR << "CSX file error\n";
//         exit(1);
//     } catch (boost::archive::archive_exception &e) {
//         LOG_ERROR << "dumping CSX to file failed: " << e.what() << "\n";
//         exit(1);
//     }
//     //file.close(); //Any open file is automatically closed when the ofstream is destroyed
// }

// template<typename IndexType, typename ValueType>
// spm_mt_t *RestoreCsx(const char *filename, IndexType **permutation)
// {
//     unsigned int nr_threads;
//     bool symmetric, reordered, full_column_indices = false;
//     ifstream file;
//     spm_mt_t *spm_mt = 0;
//     csx_sym_t<ValueType> *csx_sym = 0;
//     csx_t<ValueType> *csx = 0;

// #ifdef SPM_NUMA
//     NumaAllocator &numa_alloc = NumaAllocator::GetInstance();
// #endif

//     if (!filename) {
//         string default_filename(boost::archive::tmpdir());
//         default_filename += "/csx_file";
//         filename = default_filename.c_str();
//     }

//     try {
//         file.open(filename, ios::binary | ios::in);
//         if (!file.good()) {
//             throw ios_base::failure("");
//         }

//         boost::archive::binary_iarchive ia(file);
//         ia >> nr_threads & symmetric;

//         // Construct spm_mt
//         spm_mt = new spm_mt_t;
//         spm_mt->nr_threads = nr_threads;
//         spm_mt->symmetric = symmetric;
//         spm_mt->spm_threads = new spm_mt_thread_t[nr_threads];

// #ifdef SPM_NUMA
//         full_column_indices = true;
// #endif

//         for (unsigned int i = 0; i < nr_threads; i++) {
//             ia >> spm_mt->spm_threads[i].cpu & spm_mt->spm_threads[i].id 
//                 & spm_mt->spm_threads[i].node;
// #ifdef SPM_NUMA
//             int node = spm_mt->spm_threads[i].node;
//             if (symmetric)
//                 csx_sym = new (numa_alloc, node) csx_sym_t<ValueType>;

//             csx = new (numa_alloc, node) csx_t<ValueType>;
// #else
//             if (symmetric)
//                 csx_sym = new csx_sym_t<ValueType>;

//             csx = new csx_t<ValueType>;
// #endif
//             if (!csx) {
//                 LOG_ERROR << "CSX allocation failed\n";
//                 exit(1);
//             }

//             ia >> csx->nnz & csx->ncols & csx->nrows & csx->ctl_size
//                 & csx->row_start;
// #ifdef SPM_NUMA
//             if (symmetric)
//                 csx_sym->dvalues = new (numa_alloc, node) ValueType[csx->nrows];

//             csx->values = new (numa_alloc, node) ValueType[csx->nnz];

//             int alloc_err = 0;
//             alloc_err = check_region(csx->values, csx->nnz*sizeof(*csx->values),
//                                      node);
//             print_alloc_status("values", alloc_err);

//             csx->ctl = new (numa_alloc, node) uint8_t[csx->ctl_size];
//             alloc_err = check_region(csx->ctl, csx->ctl_size*sizeof(uint8_t),
//                                      node);
//             print_alloc_status("ctl", alloc_err);
//             csx->rows_info = new (numa_alloc, node) row_info_t[csx->nrows];
// #else
//             if (symmetric)
//                 csx_sym->dvalues = new ValueType[csx->nrows];

//             csx->values = new ValueType[csx->nnz];
//             csx->ctl = new uint8_t[csx->ctl_size];
//             csx->rows_info = new row_info_t[csx->nrows];
// #endif    
//             if (!csx->values || !csx->ctl) {
//                 LOG_ERROR << "CSX allocation failed\n";
//                 exit(1);
//             }

//             ia >> boost::serialization::make_array(csx->values, csx->nnz);
//             ia >> boost::serialization::make_array(csx->ctl, csx->ctl_size);
//             ia >> csx->id_map;
//             ia >> csx->row_jumps;
//             ia >> boost::serialization::make_array(csx->rows_info, csx->nrows);

//             if (symmetric) {
//                 ia >> boost::serialization::make_array(csx_sym->dvalues, csx->nrows);
//                 map_t *map = 0;
//                 unsigned int length;
//                 ia >> length;
// #ifdef SPM_NUMA
//                 map = new (numa_alloc, node) map_t;
//                 map->cpus = new (numa_alloc, node) unsigned int[length];
//                 map->elems_pos = new (numa_alloc, node) unsigned int[length];
// #else
//                 map = new map_t;
//                 map->cpus = new unsigned int[length];
//                 map->elems_pos = new unsigned int[length];
// #endif
//                 map->length = length;
//                 ia >> boost::serialization::make_array(map->cpus, length);
//                 ia >> boost::serialization::make_array(map->elems_pos, length);
//                 spm_mt->spm_threads[i].map = map;
//                 csx_sym->lower_matrix = csx;
//                 spm_mt->spm_threads[i].spm = csx_sym;
//             } else {
//                 spm_mt->spm_threads[i].spm = csx;
//             }

//             spm_mt->spm_threads[i].row_start = csx->row_start;
//             spm_mt->spm_threads[i].nr_rows = csx->nrows;
//         }

//         ia >> reordered;
//         if (reordered) {
//             *permutation = new IndexType[csx->ncols];
//             ia >> boost::serialization::make_array(*permutation, csx->ncols);
//         }

//         // Initialize the CSX JIT execution engine
//         CsxExecutionEngine &engine = CsxJitInit();
//         CsxJit<IndexType, ValueType> **Jits =
//             new CsxJit<IndexType, ValueType>*[nr_threads];
//         for (size_t i = 0; i < nr_threads; ++i) {
//             if (symmetric) {
//                 csx_sym = (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;
//                 csx = (csx_t<ValueType> *) csx_sym->lower_matrix;
//             } else { 
//                 csx = (csx_t<ValueType> *) spm_mt->spm_threads[i].spm;
//             }
//             bool row_jumps = csx->row_jumps != 0;
//             Jits[i] = new CsxJit<IndexType, ValueType>(csx, &engine,
//                                                        i, symmetric, 
//                                                        full_column_indices,
//                                                        row_jumps);
//             Jits[i]->GenCode(std::cout);
//         }

//         for (size_t i = 0; i < nr_threads; i++)
//             spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();

//         // Cleanup
//         for (size_t i = 0; i < nr_threads; i++)
//             delete Jits[i];
//         delete[] Jits;
//     } catch (ios_base::failure &e) {
//         LOG_ERROR << "CSX file error\n";
//         exit(1);
//     } catch (boost::archive::archive_exception &e) {
//         LOG_ERROR << "loading CSX from file failed: " << e.what() << "\n";
//         exit(1);
//     } catch (bad_alloc &e) {
//         LOG_ERROR << "loading CSX from file failed: " << e.what() << "\n";
//         exit(1);
//     }
//     //file.close(); //Any open file is automatically closed when the ifstream is destroyed

//     return spm_mt;
// }
