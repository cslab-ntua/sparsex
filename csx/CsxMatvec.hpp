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
#ifndef CSX_MATVEC_HPP
#define CSX_MATVEC_HPP

#include "SpmMt.hpp"
#include "SpmvMethod.hpp"
#include "Vector.hpp"

#ifdef __cplusplus

#include "Affinity.hpp"
#include "Barrier.hpp"

#ifdef DISABLE_POOL
#   include <boost/thread/barrier.hpp>
#endif
#include <boost/thread/thread.hpp>
#include <vector>
#endif

BEGIN_C_DECLS

void MatVecMult(spm_mt_t *spm_mt, vector_t *x, spx_scalar_t alpha, 
                vector_t *y);
void MatVecMult_sym(spm_mt_t *spm_mt, vector_t *x, spx_scalar_t alpha,
                    vector_t *y);
void MatVecKernel(spm_mt_t *spm_mt, vector_t *x, spx_scalar_t alpha, 
                  vector_t *y, spx_scalar_t beta);
void MatVecKernel_sym(spm_mt_t *spm_mt, vector_t *x, spx_scalar_t alpha,
                      vector_t *y, spx_scalar_t beta);

END_C_DECLS

void do_mv_thread(void *args);
void do_mv_sym_thread(void *args);
void do_kernel_thread(void *params);
void do_kernel_sym_thread(void *args);

#endif // CSX_MATVEC_HPP
