/*
 * CsxMatvec.hpp -- Multithreaded kernel y <-- alpha*A*x + beta*y
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2013,      Athena Elafrou
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
#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <vector>

extern "C" {
#endif

void do_matvec_thread(void *params);
void MatVecMult(spm_mt_t *spm_mt, spx_vector_t *x, spx_scalar_t alpha, 
                spx_vector_t *y, spx_scalar_t beta);
void MatVecMult_sym(spm_mt_t *spm_mt, spx_vector_t *x, spx_scalar_t alpha, 
                    spx_vector_t *y, spx_scalar_t beta);

#ifdef __cplusplus
}
#endif

#endif // CSX_MATVEC_HPP
