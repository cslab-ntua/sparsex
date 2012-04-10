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

/* Elmer types for CSR implementation */
typedef int32_t elmer_index_t;
typedef double elmer_value_t;

void elmer_matvec_(void **tuned, void *n, void *rowptr, void *colind,
                   void *values, void *x, void *y, void *reinit);

#endif  /* ELMERIF_H__ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
