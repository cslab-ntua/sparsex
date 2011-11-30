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

#ifndef DELTA_TEMPLATE
#   define DELTA_TEMPLATE MULT_TEMPLATE_DIR "/delta_tmpl.c"
#endif

#ifndef HORIZ_TEMPLATE
#   define HORIZ_TEMPLATE MULT_TEMPLATE_DIR "/horiz_tmpl.c"
#endif

const std::string CsxTemplateSource = CSX_TEMPLATE;
const std::string DeltaTemplateSource = DELTA_TEMPLATE;
const std::string HorizTemplateSource = HORIZ_TEMPLATE;

#endif // JIT_CONFIG_H__
