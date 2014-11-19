/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file PoskiModule.hpp
 * \brief The SpMV kernel implemented with pOSKI
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef POSKI_MODULE_HPP
#define POSKI_MODULE_HPP

#include <sparsex/internals/cdecl.h>

SPX_BEGIN_C_DECLS__
// Undefine HAVE_CONFIG_H, as if we were including poski normally
#undef HAVE_CONFIG_H
#include <poski/poski.h>
SPX_END_C_DECLS__

#include "Timer.hpp"
#include <sparsex/types.h>
#include <iostream>

void poski_spmv(spx_index_t *Aptr, spx_index_t *Aind, spx_value_t *Aval,
                spx_index_t nrows, spx_index_t ncols, spx_index_t nnz,
                spx_value_t *x, spx_value_t *y,
                spx_value_t ALPHA, spx_value_t BETA);

#endif  // POSKI_MODULE_HPP

