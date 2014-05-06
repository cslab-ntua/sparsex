/*
 * Map.hpp -- Map used for symmetric sparse matrices.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef MAP_HPP
#define MAP_HPP

///> The definition of map used for symmetric sparse matrices.
typedef struct map {
    unsigned int length;
    unsigned int *cpus;
    unsigned int *elems_pos;
} map_t;

#endif /* MAP_HPP */
