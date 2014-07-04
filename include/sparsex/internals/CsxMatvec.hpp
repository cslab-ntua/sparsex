/*
 * CsxMatvec.hpp -- Multithreaded kernel y <-- alpha*A*x + beta*y
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_MATVEC_HPP
#define SPARSEX_INTERNALS_CSX_MATVEC_HPP

#include "sparsex/internals/SpmMt.hpp"
#include "sparsex/internals/SpmvMethod.hpp"
#include "sparsex/internals/Vector.hpp"

#ifdef __cplusplus

#include "sparsex/internals/Affinity.hpp"
#include "sparsex/internals/Barrier.hpp"

#ifdef DISABLE_POOL
#   include <boost/thread/barrier.hpp>
#endif
#include <boost/thread/thread.hpp>
#include <vector>
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

void do_mv_thread(void *args);
void do_mv_sym_thread(void *args);
void do_kernel_thread(void *params);
void do_kernel_sym_thread(void *args);

#endif // SPARSEX_INTERNALS_CSX_MATVEC_HPP