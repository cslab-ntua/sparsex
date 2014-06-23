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

#include "sparsex/internals/CtlUtil.hpp"

///< CSX matrix format
typedef struct {
    unsigned long rowptr;    /* rowptr is the index in csx->ctl of
                                the first element of row i */
    unsigned long valptr;    /* valptr is the index in csx->values of
                                the first element of row i */
    unsigned long span;
} row_info_t;

typedef struct {
    unsigned long nnz, ncols, nrows, ctl_size, row_start;
    uint8_t row_jumps;
    double *values;
    uint8_t *ctl;
    long id_map[CTL_PATTERNS_MAX];
    row_info_t *rows_info; 
} csx_double_t;

typedef struct {
    csx_double_t *lower_matrix;
    double *dvalues;
} csx_double_sym_t;

#ifdef __cplusplus
// C++ only

template<typename ValueType>
struct csx_t {
    unsigned long nnz, ncols, nrows, ctl_size, row_start;
    uint8_t row_jumps;
    ValueType *values;
    uint8_t *ctl;
    long id_map[CTL_PATTERNS_MAX];
    row_info_t *rows_info;
};

template<typename ValueType>
struct csx_sym_t {
    csx_t<ValueType> *lower_matrix;
    ValueType *dvalues;
};

#endif

#endif  // SPARSEX_INTERNALS_CSX_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
