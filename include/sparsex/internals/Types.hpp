/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Types.hpp
 * \brief Wrap API's types.h and promote user supplied types to our internals
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_TYPES_HPP
#define SPARSEX_INTERNALS_TYPES_HPP

#include <sparsex/types.h>
#include <boost/type_traits.hpp>

// Index type should be promoted to unsigned internally
typedef boost::make_unsigned<spx_index_t>::type spx_uindex_t;

#endif  // SPARSEX_INTERNALS_TYPES_HPP
