/* -*- C++ -*-
 *
 * CsxUtil.hpp -- CSX-related routines.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_UTIL_HPP
#define CSX_UTIL_HPP

#include "Allocators.hpp"
#include "Csx.hpp"
#include "SpmMt.hpp"
#include <cstdlib>

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX form.
 *
 *  @param spm_mt  the sparse matrix in CSX format.
 */
template<typename ValueType>
uint64_t CsxSize(void *spm_mt);

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX-Sym form.
 *
 *  @param spm_mt  the sparse matrix in CSX-Sym format.
 */
template<typename ValueType>
uint64_t CsxSymSize(void *spm_mt);


/**
 *  Find the size of the map including the values accessed by it.
 *
 *  @param spm_mt  the complete multithreaded CSX-Sym matrix.
 *  @return        the size of the map.
 */
uint64_t MapSize(void *spm);

/**
 *  Deallocation of CSX or CSX-Sym sparse matrix.
 *  
 *  @param spm_mt the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
template<typename ValueType>
void PutSpmMt(spm_mt_t *spm_mt);

///> Deletes all the fields of CSX format.
template<typename ValueType>
void DestroyCsx(csx_t<ValueType> *csx);

///> Deletes all the fields of CSX-Sym format.
template<typename ValueType>
void DestroyCsxSym(csx_sym_t<ValueType> *csx_sym);


/* Function definitions */
template<typename ValueType>
uint64_t CsxSize(void *spm_mt)
{
    return (uint64_t) CsxSize<ValueType>((spm_mt_t *) spm_mt);
}

template<typename ValueType>
unsigned long CsxSize(spm_mt_t *spm_mt)
{
    unsigned long ret;

    ret = 0;
    for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        csx_t<ValueType> *csx = (csx_t<ValueType> *)t->spm;
        
        ret += csx->nnz * sizeof(ValueType);
        ret += csx->ctl_size;
    }

    return ret;
}

template<typename ValueType>
uint64_t CsxSymSize(void *spm_mt)
{
    return (uint64_t) CsxSymSize<ValueType>((spm_mt_t *) spm_mt);
}

template<typename ValueType>
unsigned long CsxSymSize(spm_mt_t *spm_mt)
{
    unsigned long ret;

    ret = 0;
    for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        csx_sym_t<ValueType> *csx_sym = (csx_sym_t<ValueType> *)t->spm;
        csx_t<ValueType> *csx = csx_sym->lower_matrix;

        ret += (csx->nnz + csx->nrows) * sizeof(ValueType);
        ret += csx->ctl_size;
    }

    return ret;
}    

template<typename ValueType>
void PutSpmMt(spm_mt_t *spm_mt)
{
#ifdef SPM_NUMA
        if (spm_mt->interleaved) {
            size_t total_ctl_size = 0, total_nnz = 0;
            csx_t<ValueType> *csx = 0;
            for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
                if (spm_mt->symmetric) {
                    csx_sym_t<ValueType> *csx_sym =
                        (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;           
                    csx = csx_sym->lower_matrix;
                } else {
                    csx = (csx_t<ValueType> *) spm_mt->spm_threads[i].spm;
                }

                total_nnz += csx->nnz;
                total_ctl_size += csx->ctl_size;
                if (i != 0) {
                    csx->ctl = 0;
                    csx->values = 0;
                }
            }

            NumaAllocator &alloc = NumaAllocator::GetInstance();
            if (spm_mt->symmetric) {
                csx_sym_t<ValueType> *csx_sym =
                    (csx_sym_t<ValueType> *) spm_mt->spm_threads[0].spm;
                csx = csx_sym->lower_matrix;
            } else {
                csx = (csx_t<ValueType> *) spm_mt->spm_threads[0].spm;
            }

            alloc.Destroy(csx->ctl, total_ctl_size);
            alloc.Destroy(csx->values, total_nnz);
        }
#endif

    if (!spm_mt->symmetric) {
        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            csx_t<ValueType> *csx = (csx_t<ValueType> *)
                spm_mt->spm_threads[i].spm;
            DestroyCsx(csx);
        }
    } else {
        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            csx_sym_t<ValueType> *csx_sym =
                (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;           
            DestroyCsxSym(csx_sym);
#ifdef SPM_NUMA
            NumaAllocator &alloc = NumaAllocator::GetInstance();
            alloc.Destroy(spm_mt->spm_threads[i].map->cpus, 
                          spm_mt->spm_threads[i].map->length);
            alloc.Destroy(spm_mt->spm_threads[i].map->elems_pos,
                          spm_mt->spm_threads[i].map->length);
            alloc.Destroy(spm_mt->spm_threads[i].map);
#else
            delete[] spm_mt->spm_threads[i].map->cpus;
            delete[] spm_mt->spm_threads[i].map->elems_pos;
            delete spm_mt->spm_threads[i].map;
#endif
        }
    }

    delete[] spm_mt->spm_threads;
    delete spm_mt;
}

template<typename ValueType>
static void DestroyCsx(csx_t<ValueType> *csx)
{
#ifdef SPM_NUMA
    NumaAllocator &alloc = NumaAllocator::GetInstance();
    if (csx->ctl && csx->values) {
        alloc.Destroy(csx->ctl, csx->ctl_size);
        alloc.Destroy(csx->values, csx->nnz);
    }
    alloc.Destroy(csx->rows_info, csx->nrows);
    alloc.Destroy(csx);
#else
    delete[] csx->ctl;
    delete[] csx->values;
    delete[] csx->rows_info;
    delete csx;
#endif
}

template<typename ValueType>
static void DestroyCsxSym(csx_sym_t<ValueType> *csx_sym)
{
#ifdef SPM_NUMA
    uint64_t diag_size = csx_sym->lower_matrix->nrows;
#endif

    DestroyCsx(csx_sym->lower_matrix);
#ifdef SPM_NUMA
    NumaAllocator &alloc = NumaAllocator::GetInstance();
    alloc.Destroy(csx_sym->dvalues, diag_size);
    alloc.Destroy(csx_sym);
#else
    delete[] csx_sym->dvalues;
    delete csx_sym;
#endif
}

#endif  // CSX_UTIL_HPP
