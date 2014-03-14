/*
 * CsxSpmvMt.hpp
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SPMV_MT_HPP
#define CSX_SPMV_MT_HPP

#include "Affinity.hpp"
#include "Csr.hpp"
#include "SpmMt.hpp"
#include "SpmvMethod.hpp"
#include "Types.hpp"
#include "Vector.hpp"
#include "tsc.h"

#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

float spmv_bench_mt(spm_mt_t *spm_mt, unsigned long loops,
                    unsigned long rows_nr, unsigned long cols_nr);
void spmv_check_mt(csx::CSR<index_t, value_t> *spm, spm_mt_t *spm_mt,
                   unsigned long loops, unsigned long rows_nr,
                   unsigned long cols_nr);
float spmv_bench_sym_mt(spm_mt_t *spm_mt, unsigned long loops,
                        unsigned long rows_nr, unsigned long cols_nr);
void spmv_check_sym_mt(csx::CSR<index_t, value_t> *spm, spm_mt_t *spm_mt,
                       unsigned long loops, unsigned long rows_nr,
                       unsigned long cols_nr);

#endif  // CSX_SPMV_MT_HPP

