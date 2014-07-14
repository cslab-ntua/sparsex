/*
 * \file CsxKernels.hpp
 *
 * \brief Multithreaded SpMV kernels
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_KERNELS_HPP
#define SPARSEX_INTERNALS_CSX_KERNELS_HPP

#include <sparsex/internals/SpmMt.hpp>
#include <sparsex/internals/Vector.hpp>

#ifdef __cplusplus

#include <sparsex/internals/Affinity.hpp>
#include <sparsex/internals/CsxSpmv.hpp>

#if SPX_DISABLE_POOL
#   include <boost/thread/barrier.hpp>
#   include <boost/thread/thread.hpp>
#   include <vector>
#else
#   include <sparsex/internals/Barrier.hpp>
#   include <sparsex/internals/ThreadPool.hpp>
#endif

#endif

SPX_BEGIN_C_DECLS__

void MatVecMult(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha, 
                vector_t *y);
void MatVecMult_sym(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha,
                    vector_t *y);
void MatVecKernel(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha, 
                  vector_t *y, spx_value_t beta);
void MatVecKernel_sym(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha,
                      vector_t *y, spx_value_t beta);

SPX_END_C_DECLS__

#endif // SPARSEX_INTERNALS_CSX_KERNELS_HPP
