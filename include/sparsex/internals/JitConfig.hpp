/*
 * JitConfig.hpp -- JIT global configuration options.
 *
 * Copyright (C) 2011-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_JIT_CONFIG_HPP
#define SPARSEX_INTERNALS_JIT_CONFIG_HPP

#include <sparsex/internals/Config.hpp>

using namespace std;

#ifndef CSX_TEMPLATE
#   define CSX_TEMPLATE SPX_MULT_TEMPLATE_DIR "/csx_spmv_tmpl.c"
#endif

#ifndef CSX_SYM_TEMPLATE
#   define CSX_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/csx_sym_spmv_tmpl.c"
#endif

#ifndef DELTA_TEMPLATE
#   define DELTA_TEMPLATE SPX_MULT_TEMPLATE_DIR "/delta_tmpl.c"
#endif

#ifndef HORIZ_TEMPLATE
#   define HORIZ_TEMPLATE SPX_MULT_TEMPLATE_DIR "/horiz_tmpl.c"
#endif

#ifndef VERT_TEMPLATE
#   define VERT_TEMPLATE SPX_MULT_TEMPLATE_DIR "/vert_tmpl.c"
#endif

#ifndef DIAG_TEMPLATE
#   define DIAG_TEMPLATE SPX_MULT_TEMPLATE_DIR "/diag_tmpl.c"
#endif

#ifndef RDIAG_TEMPLATE
#   define RDIAG_TEMPLATE SPX_MULT_TEMPLATE_DIR "/rdiag_tmpl.c"
#endif

#ifndef BLOCK_ROW_TEMPLATE
#   define BLOCK_ROW_TEMPLATE SPX_MULT_TEMPLATE_DIR "/block_row_tmpl.c"
#endif

#ifndef BLOCK_ROW_ONE_TEMPLATE
#   define BLOCK_ROW_ONE_TEMPLATE SPX_MULT_TEMPLATE_DIR "/block_row_one_tmpl.c"
#endif

#ifndef BLOCK_COL_TEMPLATE
#   define BLOCK_COL_TEMPLATE SPX_MULT_TEMPLATE_DIR "/block_col_tmpl.c"
#endif

#ifndef BLOCK_COL_ONE_TEMPLATE
#   define BLOCK_COL_ONE_TEMPLATE SPX_MULT_TEMPLATE_DIR "/block_col_one_tmpl.c"
#endif

#ifndef DELTA_SYM_TEMPLATE
#   define DELTA_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/delta_sym_tmpl.c"
#endif

#ifndef HORIZ_SYM_TEMPLATE
#   define HORIZ_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/horiz_sym_tmpl.c"
#endif

#ifndef VERT_SYM_TEMPLATE
#   define VERT_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/vert_sym_tmpl.c"
#endif

#ifndef DIAG_SYM_TEMPLATE
#   define DIAG_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/diag_sym_tmpl.c"
#endif

#ifndef RDIAG_SYM_TEMPLATE
#   define RDIAG_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/rdiag_sym_tmpl.c"
#endif

#ifndef BLOCK_ROW_SYM_TEMPLATE
#   define BLOCK_ROW_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/block_row_sym_tmpl.c"
#endif

#ifndef BLOCK_COL_SYM_TEMPLATE
#   define BLOCK_COL_SYM_TEMPLATE SPX_MULT_TEMPLATE_DIR "/block_col_sym_tmpl.c"
#endif

const string CsxTemplateSource = CSX_TEMPLATE;
const string CsxSymTemplateSource = CSX_SYM_TEMPLATE;
const string DeltaTemplateSource = DELTA_TEMPLATE;
const string HorizTemplateSource = HORIZ_TEMPLATE;
const string VertTemplateSource = VERT_TEMPLATE;
const string DiagTemplateSource = DIAG_TEMPLATE;
const string RDiagTemplateSource = RDIAG_TEMPLATE;
const string BlockRowTemplateSource = BLOCK_ROW_TEMPLATE;
const string BlockRowOneTemplateSource = BLOCK_ROW_ONE_TEMPLATE;
const string BlockColTemplateSource = BLOCK_COL_TEMPLATE;
const string BlockColOneTemplateSource = BLOCK_COL_ONE_TEMPLATE;
const string DeltaSymTemplateSource = DELTA_SYM_TEMPLATE;
const string HorizSymTemplateSource = HORIZ_SYM_TEMPLATE;
const string VertSymTemplateSource = VERT_SYM_TEMPLATE;
const string DiagSymTemplateSource = DIAG_SYM_TEMPLATE;
const string RDiagSymTemplateSource = RDIAG_SYM_TEMPLATE;
const string BlockRowSymTemplateSource = BLOCK_ROW_SYM_TEMPLATE;
const string BlockColSymTemplateSource = BLOCK_COL_SYM_TEMPLATE;

#endif // SPARSEX_INTERNALS_JIT_CONFIG_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
