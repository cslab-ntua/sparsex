/*
 * Copyright (C) 2011-2016, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2016, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file futex.h
 * \brief Linux futex wrapper
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2016
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_FUTEX_H
#define SPARSEX_INTERNALS_FUTEX_H

#include <sparsex/internals/cdecl.h>

SPX_BEGIN_C_DECLS__

void futex_wait(int *addr, int val);
void futex_wake(int *addr, int count);

SPX_END_C_DECLS__

#endif  // SPARSEX_INTERNALS_FUTEX_H
