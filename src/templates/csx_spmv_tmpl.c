/*
 * Copyright (C) 2009-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */ 

/**
 * \file csx_spmv_tmpl.c
 * \brief The CSX multiplication template for runtime code generation.
 * 
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/types.h>
#include <sparsex/internals/CtlUtil.hpp>
#include <sparsex/internals/Vector.hpp>
#include <sparsex/internals/Csx.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#define CSX_SPMV_FN_MAX CTL_PATTERNS_MAX

#define ALIGN(buf, a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
// the following function will be unused when only delta8 units
// are present
static void align_ptr(uint8_t **ctl, int align)
{
  *ctl = ALIGN(*ctl, align);
}
#pragma GCC diagnostic pop

#if SPX_DEBUG
static void ctl_print(uint8_t *ctl, spx_index_t start, spx_index_t end,
                      const char *descr)
{
  for (spx_index_t i = start; i < end; i++) {
    printf("%s[%ld]: %p = %d\n", descr, i, &ctl[i], ctl[i]);
    fflush(stdout);
  }
}

static void deref(void *ptr)
{
  volatile unsigned long val = *((unsigned long *) ptr);
  val++;
}
#endif

typedef spx_value_t (csx_spmv_fn_t)(uint8_t **ctl, uint8_t size,
                                    spx_value_t **values,
                                    spx_value_t **x, spx_value_t **y);

${spmv_func_definitions}

void spm_csx_multiply(void *spm, vector_t *in, vector_t *out,
                      spx_value_t scale_f, vector_t *local_out /* unused */)
{
  csx_matrix_t *csx = (csx_matrix_t *) spm;
  spx_value_t *x = in->elements;
  spx_value_t *y = out->elements;
  spx_value_t *v = csx->values;
  spx_value_t *x_curr = x;
  spx_value_t *y_curr = y + csx->row_start;
  register spx_value_t yr = 0;
  uint8_t *ctl = csx->ctl;
  uint8_t *ctl_end = ctl + csx->ctl_size;
  uint8_t size, flags;
  uint8_t patt_id;

  /* uint8_t *ctl_start = ctl; */
  /* ctl_print(ctl, 0, csx->ctl_size, "ctl"); */
  do {
    flags = *ctl++;
    size = *ctl++;
    if (test_bit(&flags, CTL_NR_BIT)) {
      *y_curr += yr;
      yr = 0;
      ${new_row_hook}
      x_curr = x;
    }

    ${next_x}
    patt_id = flags & CTL_PATTERN_MASK;
    ${body_hook}
    /* printf("ctl moved at %zd bytes\n", ctl - ctl_start); */
    /* fflush(stdout); */
  } while (ctl < ctl_end);

  *y_curr += yr;
}
