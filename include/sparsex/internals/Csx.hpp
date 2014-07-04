/*
 * Csx.hpp -- The CSX data structure (C front-end)
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

#ifndef SPARSEX_INTERNALS_CSX_HPP
#define SPARSEX_INTERNALS_CSX_HPP

#include <sparsex/types.h>
#include <sparsex/internals/CtlUtil.hpp>

///< CSX matrix format
typedef struct {
    spx_index_t rowptr;    /* rowptr is the index in csx->ctl of
                                the first element of row i */
    spx_index_t valptr;    /* valptr is the index in csx->values of
                                the first element of row i */
    spx_index_t span;
} row_info_t;

typedef struct {
    spx_value_t *values;
    spx_index_t nnz;
    spx_index_t ncols;
    spx_index_t nrows;
    spx_index_t ctl_size;
    spx_index_t row_start;
    uint8_t *ctl;
    uint8_t row_jumps;
    long id_map[CTL_PATTERNS_MAX];
    row_info_t *rows_info; 
} csx_matrix_t;

typedef struct {
    csx_matrix_t *lower_matrix;
    spx_value_t *dvalues;
} csx_sym_matrix_t;

#ifdef __cplusplus
// C++ only

// FIXME: members of CsxMatrix MUST have the same order with csx_matrix_t
//        This is error-prone with explicit casts between the two types
template<typename IndexType, typename ValueType>
struct CsxMatrix {
    ValueType *values;
    IndexType nnz;
    IndexType ncols;
    IndexType nrows;
    IndexType ctl_size;
    IndexType row_start;
    uint8_t *ctl;
    uint8_t row_jumps;
    long id_map[CTL_PATTERNS_MAX];
    row_info_t *rows_info;
};

template<typename IndexType, typename ValueType>
struct CsxSymMatrix {
    CsxMatrix<IndexType, ValueType> *lower_matrix;
    ValueType *dvalues;
};

#endif

#endif  // SPARSEX_INTERNALS_CSX_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
