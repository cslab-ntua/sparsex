/*
 * \file CsxUtil.hpp
 *
 * \breif CSX-related utilities
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

#ifndef SPARSEX_INTERNALS_CSX_UTIL_HPP
#define SPARSEX_INTERNALS_CSX_UTIL_HPP

#include <sparsex/internals/Allocators.hpp>
#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Csx.hpp>
#include <sparsex/internals/Element.hpp>
#include <sparsex/internals/SpmMt.hpp>
#include <sparsex/internals/Vector.hpp>
#include <cstdlib>

namespace sparsex {
namespace csx {

// Pattern ID generation
extern const unsigned long DeltaIdOffset;
///< ID offset for delta units.

extern const unsigned long PatternIdOffset;
///< ID offset for substructure units.

/**
 *  Generate pattern id for delta units
 *
 *  @param delta_size   byte count of delta unit
 *  @return the pattern id of the delta unit
 */
unsigned long GetPatternId(size_t delta_size);

/**
 *  Generate pattern id for CSX units
 *
 *  @param elem a CSX generic element; must be a pattern
 *  @return the pattern id of the CSX unit
 *  @throws invalid_argument if elem is not a pattern
 */
template<typename IndexType, typename ValueType>
unsigned long GetPatternId(const Element<IndexType, ValueType> &elem)
{
    const Encoding::Instantiation &inst = elem.GetInstantiation();
    size_t size = elem.GetSize();
    if (inst.first == Encoding::None)
        throw invalid_argument("elem is not a pattern");

    Encoding e(inst.first);
    unsigned long ret;
    if (e.IsBlock())
        ret = inst.first*PatternIdOffset + size / e.GetBlockAlignment();
    else
        ret = inst.first*PatternIdOffset + inst.second;

    return ret;
}

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
 *  Find the size of the map including the values accessed by it.
 *
 *  @param spm_mt  the complete multithreaded CSX-Sym matrix.
 *  @return        the size of the map.
 */
template<typename IndexType, typename ValueType>
size_t CsxSymMapSize(void *spm);

/**
 *  Deallocation of CSX or CSX-Sym sparse matrix.
 *  
 *  @param spm_mt the (multithreaded) CSX or CSX-Sym sparse matrix.
 */
template<typename IndexType, typename ValueType>
void PutSpmMt(spm_mt_t *spm_mt);

///> Deletes all the fields of CSX format.
template<typename IndexType, typename ValueType>
void DestroyCsx(CsxMatrix<IndexType, ValueType> *csx);

///> Deletes all the fields of CSX-Sym format.
template<typename IndexType, typename ValueType>
void DestroyCsxSym(CsxSymMatrix<IndexType, ValueType> *csx_sym);


/* Function definitions */
template<typename IndexType, typename ValueType>
size_t CsxSize(void *spm_mt)
{
    return CsxSize<IndexType, ValueType>((spm_mt_t *) spm_mt);
}

template<typename IndexType, typename ValueType>
size_t CsxSize(spm_mt_t *spm_mt)
{
    size_t ret = 0;
    for (size_t i = 0; i < spm_mt->nr_threads; ++i) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        CsxMatrix<IndexType, ValueType> *csx =
            (CsxMatrix<IndexType, ValueType> *) t->spm;
        
        ret += csx->nnz * sizeof(ValueType);
        ret += csx->ctl_size;
    }

    return ret;
}

template<typename IndexType, typename ValueType>
size_t CsxSymSize(void *spm_mt)
{
    return CsxSymSize<IndexType, ValueType>((spm_mt_t *) spm_mt);
}

template<typename IndexType, typename ValueType>
size_t CsxSymSize(spm_mt_t *spm_mt)
{
    size_t ret;

    ret = 0;
    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        CsxSymMatrix<IndexType, ValueType> *csx_sym =
            (CsxSymMatrix<IndexType, ValueType> *)t->spm;
        CsxMatrix<IndexType, ValueType> *csx = csx_sym->lower_matrix;

        ret += (csx->nnz + csx->nrows) * sizeof(ValueType);
        ret += csx->ctl_size;
    }

    return ret;
}    

template<typename IndexType, typename ValueType>
size_t CsxSymMapSize(void *spm)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread;
    size_t size = 0;
    
    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        size += spm_thread->map->length * sizeof(IndexType);
        size += spm_thread->map->length * sizeof(IndexType);
        size += spm_thread->map->length * sizeof(ValueType);
    }
    
    return size;
}

template<typename IndexType, typename ValueType>
void PutSpmMt(spm_mt_t *spm_mt)
{
#if SPX_USE_NUMA
        if (spm_mt->interleaved) {
            size_t total_ctl_size = 0, total_nnz = 0;
            CsxMatrix<IndexType, ValueType> *csx = 0;
            for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
                if (spm_mt->symmetric) {
                    CsxSymMatrix<IndexType, ValueType> *csx_sym =
                        (CsxSymMatrix<IndexType, ValueType> *)
                        spm_mt->spm_threads[i].spm;           
                    csx = csx_sym->lower_matrix;
                } else {
                    csx = (CsxMatrix<IndexType, ValueType> *)
                        spm_mt->spm_threads[i].spm;
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
                CsxSymMatrix<IndexType, ValueType> *csx_sym =
                    (CsxSymMatrix<IndexType, ValueType> *)
                    spm_mt->spm_threads[0].spm;
                csx = csx_sym->lower_matrix;
            } else {
                csx = (CsxMatrix<IndexType, ValueType> *)
                    spm_mt->spm_threads[0].spm;
            }

            alloc.Destroy(csx->ctl, total_ctl_size);
            alloc.Destroy(csx->values, total_nnz);
        }
#endif

    if (!spm_mt->symmetric) {
        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            CsxMatrix<IndexType, ValueType> *csx =
                (CsxMatrix<IndexType, ValueType> *) spm_mt->spm_threads[i].spm;
            DestroyCsx(csx);
        }
    } else {
        for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
            CsxSymMatrix<IndexType, ValueType> *csx_sym =
                (CsxSymMatrix<IndexType, ValueType> *) spm_mt->spm_threads[i].spm;
            DestroyCsxSym(csx_sym);
#if SPX_USE_NUMA
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
            if (spm_mt->local_buffers && i != 0) 
                VecDestroy(spm_mt->local_buffers[i]);
        }

        free(spm_mt->local_buffers);
    }

    delete[] spm_mt->spm_threads;
    delete spm_mt;
}

template<typename IndexType, typename ValueType>
static void DestroyCsx(CsxMatrix<IndexType, ValueType> *csx)
{
#if SPX_USE_NUMA
    NumaAllocator &alloc = NumaAllocator::GetInstance();
    if (csx->ctl && csx->values) {
        alloc.Destroy(csx->ctl, csx->ctl_size);
        alloc.Destroy(csx->values, csx->nnz);
    }
    alloc.Destroy(csx->rows_info, csx->nrows);
    alloc.Destroy(csx);
#else
    StdAllocator &alloc = StdAllocator::GetInstance();
    alloc.Destroy(csx->ctl, csx->ctl_size);
    delete[] csx->values;
    delete[] csx->rows_info;
    delete csx;
#endif
}

template<typename IndexType, typename ValueType>
static void DestroyCsxSym(CsxSymMatrix<IndexType, ValueType> *csx_sym)
{
#if SPX_USE_NUMA
    uint64_t diag_size = csx_sym->lower_matrix->nrows;
#endif

    DestroyCsx(csx_sym->lower_matrix);
#if SPX_USE_NUMA
    NumaAllocator &alloc = NumaAllocator::GetInstance();
    alloc.Destroy(csx_sym->dvalues, diag_size);
    alloc.Destroy(csx_sym);
#else
    delete[] csx_sym->dvalues;
    delete csx_sym;
#endif
}

} // end of namespace csx
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_CSX_UTIL_HPP
