/*
 * export.h -- BCSR symbol accesible outside bcsr 'module'
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __BCSR_EXPORT_H__
#define __BCSR_EXPORT_H__

typedef struct bcsr_metadata_s {
    uint64_t br, bc;
} bcsr_metadata_t;

#endif
