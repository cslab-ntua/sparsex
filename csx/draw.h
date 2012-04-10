/*
 * draw.h -- Definitions for the drawing utility.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SPM_DRAW__
#define CSX_SPM_DRAW__

#include "spm.h"

void Draw(csx::SPM &spm, const char *filename, int row_start=0, int row_end=0,
          const int width=600, const int height=600);

#endif /* CSX_SPM_DRAW__ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
