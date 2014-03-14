/*
 * Types.hpp -- Wrap API's types.h and
 *              promote user supplied types to our internals
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef TYPES_HPP
#define TYPES_HPP

#include "types.h"

#include <boost/type_traits.hpp>

// Index type should be promoted to unsigned internally
typedef boost::make_unsigned<index_t>::type uindex_t;

#endif  // TYPES_HPP
