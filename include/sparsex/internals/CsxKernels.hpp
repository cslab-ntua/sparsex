/*
 * Copyright (C) 2011-2014, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CsxKernels.hpp
 * \brief Multithreaded SpMV kernels
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_KERNELS_HPP
#define SPARSEX_INTERNALS_CSX_KERNELS_HPP

#include <sparsex/internals/SpmMt.hpp>
#include <sparsex/internals/Vector.hpp>

#ifdef __cplusplus

#include <sparsex/internals/Affinity.hpp>
#include <sparsex/internals/CsxSpmv.hpp>
#include <sparsex/internals/Barrier.hpp>
#include <sparsex/internals/ThreadPool.hpp>

#endif

SPX_BEGIN_C_DECLS__

void MatVecMult(const spm_mt_t *spm_mt, const vector_t *x,
                spx_value_t alpha, vector_t *y);
void MatVecMult_sym(const spm_mt_t *spm_mt, const vector_t *x,
                    spx_value_t alpha, vector_t *y);
void MatVecKernel(const spm_mt_t *spm_mt, const vector_t *x,
                  spx_value_t alpha, vector_t *y, spx_value_t beta);
void MatVecKernel_sym(const spm_mt_t *spm_mt, const vector_t *x,
                      spx_value_t alpha, vector_t *y, spx_value_t beta);

SPX_END_C_DECLS__

#endif // SPARSEX_INTERNALS_CSX_KERNELS_HPP
