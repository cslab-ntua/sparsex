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

#include "../C-API/mattype.h"
#include "Csx.hpp"
#include "SpmMt.hpp"

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
    if (!spm_mt->symmetric) {
        for (unsigned i = 0; i < spm_mt->nr_threads; i++) {
            csx_t<ValueType> *csx = (csx_t<ValueType> *)
                spm_mt->spm_threads[i].spm;
            DestroyCsx(csx);
        }
    } else {
        for (unsigned i = 0; i < spm_mt->nr_threads; i++) {
            csx_sym_t<ValueType> *csx_sym =
                (csx_sym_t<ValueType> *) spm_mt->spm_threads[i].spm;           
            DestroyCsxSym(csx_sym);
//            free(spm_mt->spm_threads[i].map->cpus);
//            free(spm_mt->spm_threads[i].map->elems_pos);
//            free(spm_mt->spm_threads[i].map);
        }
    }
    free(spm_mt->spm_threads);
    free(spm_mt);
}

template<typename ValueType>
static void DestroyCsx(csx_t<ValueType> *csx)
{
#ifdef SPM_NUMA
    free_interleaved(csx->ctl, csx->ctl_size*sizeof(*csx->ctl));
    free_interleaved(csx->values, csx->nnz*sizeof(*csx->values));
    free_interleaved(csx->rows_info, csx->nrows*sizeof(*csx->rows_info));
    free_interleaved(csx, sizeof(*csx));
#else
    free(csx->ctl);
    free(csx->values);
    free(csx->rows_info);
    free(csx);
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
    numa_free(csx_sym->dvalues, diag_size*sizeof(*csx_sym->dvalues));
    numa_free(csx_sym, sizeof(*csx_sym));
#else
    free(csx_sym->dvalues);
    free(csx_sym);
#endif
}

/*
 * Explicit instantiation declarations: prevent implicit instantiations.
 * Code that would otherwise cause an implicit instantiation has to
 * use the explicit instatiation definition provided in the .cc
 */
//extern template void PutSpmMt<value_t>(spm_mt_t *);

#endif  // CSX_UTIL_HPP
