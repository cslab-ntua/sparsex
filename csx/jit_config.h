/* -*- C++ -*-
 *
 * jit_config.h -- JIT global configuration options.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef JIT_CONFIG_H__
#define JIT_CONFIG_H__

#ifndef MULT_TEMPLATE_DIR
#   define MULT_TEMPLATE_DIR "mult_templates"
#endif

#ifndef CSX_TEMPLATE
#   define CSX_TEMPLATE MULT_TEMPLATE_DIR "/csx_spmv_tmpl.c"
#endif

#ifndef CSX_SYM_TEMPLATE
#   define CSX_SYM_TEMPLATE MULT_TEMPLATE_DIR "/csx_sym_spmv_tmpl.c"
#endif

#ifndef DELTA_TEMPLATE
#   define DELTA_TEMPLATE MULT_TEMPLATE_DIR "/delta_tmpl.c"
#endif

#ifndef HORIZ_TEMPLATE
#   define HORIZ_TEMPLATE MULT_TEMPLATE_DIR "/horiz_tmpl.c"
#endif

#ifndef VERT_TEMPLATE
#   define VERT_TEMPLATE MULT_TEMPLATE_DIR "/vert_tmpl.c"
#endif

#ifndef DIAG_TEMPLATE
#   define DIAG_TEMPLATE MULT_TEMPLATE_DIR "/diag_tmpl.c"
#endif

#ifndef RDIAG_TEMPLATE
#   define RDIAG_TEMPLATE MULT_TEMPLATE_DIR "/rdiag_tmpl.c"
#endif

#ifndef BLOCK_ROW_TEMPLATE
#   define BLOCK_ROW_TEMPLATE MULT_TEMPLATE_DIR "/block_row_tmpl.c"
#endif

#ifndef BLOCK_ROW_ONE_TEMPLATE
#   define BLOCK_ROW_ONE_TEMPLATE MULT_TEMPLATE_DIR "/block_row_one_tmpl.c"
#endif

#ifndef BLOCK_COL_TEMPLATE
#   define BLOCK_COL_TEMPLATE MULT_TEMPLATE_DIR "/block_col_tmpl.c"
#endif

#ifndef BLOCK_COL_ONE_TEMPLATE
#   define BLOCK_COL_ONE_TEMPLATE MULT_TEMPLATE_DIR "/block_col_one_tmpl.c"
#endif

#ifndef DELTA_SYM_TEMPLATE
#   define DELTA_SYM_TEMPLATE MULT_TEMPLATE_DIR "/delta_sym_tmpl.c"
#endif

#ifndef HORIZ_SYM_TEMPLATE
#   define HORIZ_SYM_TEMPLATE MULT_TEMPLATE_DIR "/horiz_sym_tmpl.c"
#endif

#ifndef VERT_SYM_TEMPLATE
#   define VERT_SYM_TEMPLATE MULT_TEMPLATE_DIR "/vert_sym_tmpl.c"
#endif

#ifndef DIAG_SYM_TEMPLATE
#   define DIAG_SYM_TEMPLATE MULT_TEMPLATE_DIR "/diag_sym_tmpl.c"
#endif

#ifndef RDIAG_SYM_TEMPLATE
#   define RDIAG_SYM_TEMPLATE MULT_TEMPLATE_DIR "/rdiag_sym_tmpl.c"
#endif

#ifndef BLOCK_ROW_SYM_TEMPLATE
#   define BLOCK_ROW_SYM_TEMPLATE MULT_TEMPLATE_DIR "/block_row_sym_tmpl.c"
#endif

#ifndef BLOCK_COL_SYM_TEMPLATE
#   define BLOCK_COL_SYM_TEMPLATE MULT_TEMPLATE_DIR "/block_col_sym_tmpl.c"
#endif

const std::string CsxTemplateSource = CSX_TEMPLATE;
const std::string CsxSymTemplateSource = CSX_SYM_TEMPLATE;
const std::string DeltaTemplateSource = DELTA_TEMPLATE;
const std::string HorizTemplateSource = HORIZ_TEMPLATE;
const std::string VertTemplateSource = VERT_TEMPLATE;
const std::string DiagTemplateSource = DIAG_TEMPLATE;
const std::string RDiagTemplateSource = RDIAG_TEMPLATE;
const std::string BlockRowTemplateSource = BLOCK_ROW_TEMPLATE;
const std::string BlockRowOneTemplateSource = BLOCK_ROW_ONE_TEMPLATE;
const std::string BlockColTemplateSource = BLOCK_COL_TEMPLATE;
const std::string BlockColOneTemplateSource = BLOCK_COL_ONE_TEMPLATE;
const std::string DeltaSymTemplateSource = DELTA_SYM_TEMPLATE;
const std::string HorizSymTemplateSource = HORIZ_SYM_TEMPLATE;
const std::string VertSymTemplateSource = VERT_SYM_TEMPLATE;
const std::string DiagSymTemplateSource = DIAG_SYM_TEMPLATE;
const std::string RDiagSymTemplateSource = RDIAG_SYM_TEMPLATE;
const std::string BlockRowSymTemplateSource = BLOCK_ROW_SYM_TEMPLATE;
const std::string BlockColSymTemplateSource = BLOCK_COL_SYM_TEMPLATE;

#endif // JIT_CONFIG_H__
