/* -*- C++ -*-
 *
 * bench.h --  Benchmarking interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef BENCH_H__
#define BENCH_H__

#include "bench_util.h"

using namespace bench;

/* 
 *  Browse all the matrices available in the specified directory and solve
 *  the SpMV kernel with the specified #library or with all available libraries
 *  (LIBCSX, Intel MKL, pOSKI) if #library is not specified. The results of each
 *  solve are printed in the standard output and optionally in the provided
 *  #stats_file.
 *
 *  @param[in] directory    the directory of the MMF files. 
 *  @param[in] library      the library to be used for computing the SpMV kernel.
 *  @param[in] stats_file   an optional file for outputting benchmark results.
 *  @param[in] loops        the numbers of SpMV iterations.
 *  @param[in] outer_loops  the number of repeats of the benchmark.
 */
void Bench_Directory(const char *directory, const char *library,
                     const char *stats_file);

/* 
 *  Solve the SpMV kernel with the specified #library or with all available
 *  libraries (LIBCSX, Intel MKL, pOSKI) if #library is not specified for the
 *  given matrix in the MMF format. The results of each solve are printed in
 *  the standard output and optionally in the provided #stats_file.
 *
 *  @param[in] mmf_file     the MMF file where the matrix is stored. 
 *  @param[in] library      the library to be used for computing the SpMV kernel.
 *  @param[in] stats_file   an optional file for outputting benchmark results.
 *  @param[in] loops        the numbers of SpMV iterations.
 *  @param[in] outer_loops  the number of repeats of the benchmark.
 */
void Bench_Matrix(const char *mmf_file, const char *library,
                  const char *stats_file);

void Bench_Matrix(const char *mmf_file, SpmvFn fn, const char *stats_file);

#endif  // BENCH_H__
