/*
 * Copyright (C) 2012-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CsxSpmv.hpp
 * \brief Thread-level SpMV routines
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_SPMV_HPP
#define SPARSEX_INTERNALS_CSX_SPMV_HPP

#include <sparsex/internals/Affinity.hpp>
#include <sparsex/internals/Barrier.hpp>
#include <sparsex/internals/SpmMt.hpp>
#include <sparsex/internals/SpmvMethod.hpp>
#include <sparsex/internals/Vector.hpp>

namespace sparsex {

void do_mv_thread(void *args);
void do_mv_sym_thread(void *args);
void do_kernel_thread(void *params);
void do_kernel_sym_thread(void *args);

} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_CSX_SPMV_HPP
