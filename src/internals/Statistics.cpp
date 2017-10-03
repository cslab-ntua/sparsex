/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Statistics.cpp
 * \brief SparseX statistics collection framework
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/Statistics.hpp>

using namespace std;

namespace sparsex {
  namespace csx {

    typedef typename StatsCollection<StatsData>::InstStats InstStats;

    void SplitBlockData(size_t fixed_dim, size_t var_dim, size_t max_var_dim,
			const StatsData &data, InstStats &stats)
    {
      size_t nr_chunks = var_dim / max_var_dim;
      size_t rem_dim = var_dim % max_var_dim;
      size_t max_block_size = max_var_dim * fixed_dim;
      size_t nr_max_blocks = nr_chunks*data.GetPatterns();
      size_t rem_nr_nzeros = data.GetEncoded() - nr_max_blocks*max_block_size;
      stats[max_var_dim] += StatsData(nr_max_blocks*max_block_size,
				      nr_max_blocks, 0);
      if (rem_dim >= 2)
        // ignore one-dimensional remainders
        stats[rem_dim] += StatsData(rem_nr_nzeros, data.GetPatterns(), 0);
    }

    void EraseDims(InstStats &stats, vector<size_t> &dims)
    {
      for (size_t i = 0; i < dims.size(); ++i)
        stats.erase(dims[i]);
      dims.clear();
    }
    

    bool BlockSplitter::Manipulate(Encoding::Type type, InstStats &stats)
    {
      Encoding e(type);
      if (!e.IsBlock())
        return false;

      size_t fixed_dim = e.GetBlockAlignment();
      size_t max_block_dim = max_patt_size_ / fixed_dim;
      int ret = 0;
      vector<size_t> dims_to_erase;

      // Split stats to max_patt_size_ chunks
      InstStats::reverse_iterator i = stats.rbegin();
      for (; i != stats.rend() && i->first*fixed_dim > max_patt_size_; ++i) {
        SplitBlockData(fixed_dim, i->first, max_block_dim, i->second, stats);
        dims_to_erase.push_back(i->first);
        ++ret;
      }

      EraseDims(stats, dims_to_erase);
      InstStats::reverse_iterator j = stats.rbegin();
      for (i = stats.rbegin(); i != stats.rend(); ++i) {
        if ((i->second.GetEncoded() / (double) nr_nzeros_) < min_coverage_)
	  continue;
        for (; j != stats.rend() && j->first >= i->first &&
	       j->second.GetEncoded() / (double) nr_nzeros_ < min_coverage_;
             ++j) {
	  SplitBlockData(fixed_dim, j->first, i->first, j->second, stats);
	  dims_to_erase.push_back(j->first);
	  ++ret;
        }

      }

      EraseDims(stats, dims_to_erase);
      return ret;
    }

  } // end of namespace csx
} // end of namespace sparsex
