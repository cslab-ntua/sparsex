/*
 * Copyright (C) 2009-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2014, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CsxUtil.cpp
 * \brief CSX-related utilities
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/CsxUtil.hpp>

namespace sparsex {
  namespace csx {

    const unsigned long DeltaIdOffset = 0;
    const unsigned long PatternIdOffset = 10000;

    unsigned long GetPatternId(size_t delta_size)
    {
      return (delta_size << 3) + DeltaIdOffset;
    }

  } // end of namespace csx
} // end of namespace sparsex
