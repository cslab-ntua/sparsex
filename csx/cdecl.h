/*
 * cdecl.h -- C-related declarations
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CDECL_H
#define CDECL_H

#ifdef __cplusplus
#   define BEGIN_C_DECLS extern "C" {
#   define END_C_DECLS }
#else
#   define BEGIN_C_DECLS    /* nothing */
#   define END_C_DECLS      /* nothing */
#endif

#endif  // CDECL_H

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
