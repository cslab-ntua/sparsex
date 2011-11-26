/*
 * jit.cc -- Just In Time compilation routines.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <iostream>
#include <sstream>
#include <cassert>

#include "llvm/Analysis/Verifier.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/Target/TargetSelect.h"

#include "spm.h"
#include "csx.h"
#include "drle.h"
#include "newjit.h"
#include "jit_util.h"
#include "template_text.h"

using namespace llvm;
using namespace csx;

#ifndef CSX_TEMPLATE
#   define CSX_TEMPLATE "csx_spmv_tmpl.c"
#endif

const std::string CsxTemplateSource = CSX_TEMPLATE;

void CsxJitInit(void)
{
    InitializeNativeTarget();
}

CsxJit::CsxJit(CsxManager *csxmg, unsigned int id)
    : csxmg_(csxmg), module_(0)
{
    context_ = new LLVMContext();
    compiler_ = new ClangCompiler();
};


void CsxJit::DoNewRowHook(std::map<std::string, std::string> &hooks,
                          std::ostream &log) const
{
    hooks["new_row_hook"] = "printf(\"That's a new row.\\n\");";
}

void CsxJit::DoSpmvFnHook(std::map<std::string, std::string> &hooks,
                          std::ostream &log) const
{
    hooks["spmv_func_definitions"] =
        "void csx_delta8(uint8_t *ctl, uint8_t size, double *values, "
        "double *x, double *y) {\n"
        "\tprintf(\"csx delta8 mulitplication\\n\");\n}";
    hooks["spmv_func_entries"] = "\tcsx_delta8,";
}

void CsxJit::GenCode(std::ostream &log)
{
    // Load the template source
    TemplateText source_tmpl(SourceFromFile(CsxTemplateSource));

    // Fill in the hooks
    std::map<std::string, std::string> hooks;
    DoNewRowHook(hooks, log);
    DoSpmvFnHook(hooks, log);

    // Substitute and compile into an LLVM module
    module_ = DoCompile(source_tmpl.Substitute(hooks));

    // Optimize the code
    DoOptimize(module_);
}

Module *CsxJit::DoCompile(const std::string &source) const
{
    return compiler_->Compile(source, context_);
}

void CsxJit::DoOptimize(Module *module) const
{
    return;
}

void *CsxJit::GetSpmvFn() const
{
    return NULL;
}

int main(int argc, char **argv)
{
    CsxManager csxmg(NULL);
    CsxJit jit(&csxmg);
    jit.GenCode(std::cout);
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
