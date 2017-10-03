/*
 * Copyright (C) 2007-2013, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file SpmMt.hpp
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_SPM_MT_HPP
#define SPARSEX_INTERNALS_SPM_MT_HPP

#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Map.hpp>
#include <sparsex/internals/Vector.hpp>
#include <sparsex/internals/SpmvMethod.hpp>
#include <stdbool.h>
#include <stdint.h>

struct spm_mt_thread {
  spmv_fn_t spmv_fn;
  void *csx;
  void *data;
  map_t *map;
  // uint64_t *col_map;
  const vector_t *x;
  vector_t *y;
  int sense;
  double alpha;
  double beta;
  double secs;
  unsigned int cpu;
  unsigned int id;
  unsigned int node;
  unsigned int row_start;
  unsigned int end;
  unsigned int nr_rows;
  size_t size_assigned;
};
typedef struct spm_mt_thread spm_mt_thread_t;

struct spm_mt {
  spm_mt_thread_t *spm_threads;
  vector_t **local_buffers;
  unsigned int nr_threads;
  bool symmetric;
#if SPX_USE_NUMA
  bool interleaved;
#endif
};

typedef struct spm_mt spm_mt_t;

#endif // SPARSEX_INTERNALS_SPM_MT_HPP
