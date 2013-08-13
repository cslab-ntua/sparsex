/* -*- C++ -*-
 *
 * CsxUtil.h -- CSX-related routines.
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
#ifndef CSXUTIL_H__
#define CSXUTIL_H__

#include "csx.h"

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX form.
 *
 *  @param spm_mt  the sparse matrix in CSX format.
 */
template<typename IndexType, typename ValueType>
uint64_t CsxSize(void *spm_mt);

/**
 *  Compute the size (in bytes) of the compressed matrix in CSX-Sym form.
 *
 *  @param spm_mt  the sparse matrix in CSX-Sym format.
 */
template<typename IndexType, typename ValueType>
uint64_t CsxSymSize(void *spm_mt);

/**
 *  Deallocation of CSX or CSX-Sym sparse matrix.
 *  
 *  @param spm_mt the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
template<typename IndexType, typename ValueType>
void PutSpmMt(spm_mt_t *spm_mt);

///> Deletes all the fields of CSX format.
template<typename IndexType, typename ValueType>
static void DestroyCsx(csx_t<IndexType, ValueType> *csx);

///> Deletes all the fields of CSX-Sym format.
template<typename IndexType, typename ValueType>
static void DestroyCsxSym(csx_sym_t<IndexType, ValueType> *csx_sym);


/* Function definitions */
template<typename IndexType, typename ValueType>
uint64_t CsxSize(void *spm_mt)
{
    return (uint64_t) CsxSize<IndexType, ValueType>((spm_mt_t *) spm_mt);
}

template<typename IndexType, typename ValueType>
static unsigned long CsxSize(spm_mt_t *spm_mt)
{
    unsigned long ret;

    ret = 0;
    for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        csx_t<IndexType, ValueType> *csx = (csx_t<IndexType, ValueType> *)t->spm;
        
        ret += csx->nnz * sizeof(ValueType);
        ret += csx->ctl_size;
    }

    return ret;
}

template<typename IndexType, typename ValueType>
uint64_t CsxSymSize(void *spm_mt)
{
    return (uint64_t) CsxSymSize<IndexType, ValueType>((spm_mt_t *) spm_mt);
}

template<typename IndexType, typename ValueType>
static unsigned long CsxSymSize(spm_mt_t *spm_mt)
{
    unsigned long ret;

    ret = 0;
    for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        csx_sym_t<IndexType, ValueType> *csx_sym = 
            (csx_sym_t<IndexType, ValueType> *)t->spm;
        csx_t<IndexType, ValueType> *csx = csx_sym->lower_matrix;

        ret += (csx->nnz + csx->nrows) * sizeof(ValueType);
        ret += csx->ctl_size;
    }

    return ret;
}    

template<typename IndexType, typename ValueType>
void PutSpmMt(spm_mt_t *spm_mt)
{
    if (!spm_mt->symmetric) {
        for (unsigned i = 0; i < spm_mt->nr_threads; i++) {
            csx_t<IndexType, ValueType> *csx = (csx_t<IndexType, ValueType> *)
                spm_mt->spm_threads[i].spm;
            
            DestroyCsx(csx);
        }
    } else {
        for (unsigned i = 0; i < spm_mt->nr_threads; i++) {
            csx_sym_t<IndexType, ValueType> *csx_sym =
                (csx_sym_t<IndexType, ValueType> *) spm_mt->spm_threads[i].spm;
            
            DestroyCsxSym(csx_sym);
        }
    }
    free(spm_mt->spm_threads);
    free(spm_mt);
}

template<typename IndexType, typename ValueType>
static void DestroyCsx(csx_t<IndexType, ValueType> *csx)
{
#ifdef SPM_NUMA
    free_interleaved(csx->ctl, csx->ctl_size*sizeof(*csx->ctl));
    free_interleaved(csx->values, csx->nnz*sizeof(*csx->values));
    free_interleaved(csx, sizeof(*csx));
#else
    free(csx->ctl);
    free(csx->values);
    free(csx);
#endif
}

template<typename IndexType, typename ValueType>
static void DestroyCsxSym(csx_sym_t<IndexType, ValueType> *csx_sym)
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

#endif  // CSXUTIL_H__
