/*
 * elmerif.h -- Elmer interface
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef ELMERIF_H__
#define ELMERIF_H__

#include <stdint.h>

/* Elmer types for CSR implementation */
typedef int32_t elmer_index_t;
typedef double elmer_value_t;

void elmer_matvec_(void **tuned, void *n, void *rowptr, void *colind,
                   void *values, void *x, void *y, void *reinit);

void *csx_mattune(elmer_index_t *rowptr, elmer_index_t *colind,
                  elmer_value_t *values, elmer_index_t nr_rows, elmer_index_t nr_cols);

void csx_matvec(void *spm, elmer_value_t *x, elmer_index_t nr_x,
                elmer_value_t *y, elmer_index_t nr_y);


#endif  /* ELMERIF_H__ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
