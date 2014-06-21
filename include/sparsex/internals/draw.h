/*
 * draw.h -- Definitions for the drawing utility.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSEX_INTERNALS_DRAW_H
#define SPARSEX_INTERNALS_DRAW_H

#include "SparsePartition.h"

void Draw(csx::SparsePartition &spm, const char *filename, int row_start=0, int row_end=0,
          const int width=600, const int height=600);

#endif /* SPARSEX_INTERNALS_DRAW_H */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
