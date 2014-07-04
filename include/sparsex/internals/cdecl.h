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
#   define SPX_BEGIN_C_DECLS__    extern "C" {
#   define SPX_END_C_DECLS__      }
#else
#   define SPX_BEGIN_C_DECLS__    /* nothing */
#   define SPX_END_C_DECLS__      /* nothing */
#endif  // __cplusplus

#endif  // CDECL_H

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
