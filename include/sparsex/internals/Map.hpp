/*
 * Map.hpp -- Map used for symmetric sparse matrices.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_MAP_HPP
#define SPARSEX_INTERNALS_MAP_HPP

///> The definition of map used for symmetric sparse matrices.
typedef struct map {
    unsigned int length;
    unsigned int *cpus;
    unsigned int *elems_pos;
} map_t;

#endif /* SPARSEX_INTERNALS_MAP_HPP */
