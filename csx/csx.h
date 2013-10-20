/*
 * csx.h -- The CSX data structure (C front-end)
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_H__
#define CSX_H__

///< CSX matrix format
typedef struct {
    size_t rowptr;    /* rowptr is the index in csx->ctl of
                         the first element of row i */
    size_t valptr;    /* valptr is the index in csx->values of
                         the first element of row i */
    size_t span;
} row_info_t;

//padding issue with size_t???
typedef struct {
    size_t nnz, ncols, nrows, ctl_size, row_start;
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


#endif  // CSX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
