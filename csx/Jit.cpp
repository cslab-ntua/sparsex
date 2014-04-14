/*
 * Jit.cpp -- Just In Time compilation routines.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * Copyright (C) 2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "Jit.hpp"

using namespace llvm;
using namespace csx;

CsxExecutionEngine &csx::CsxJitInit(void)
{
    InitializeNativeTarget();
    return CsxExecutionEngine::CreateEngine();
}

void CsxExecutionEngine::AddModule(Module *mod)
{
    if (!llvm_engine_) {
        std::string errmsg;
        llvm_engine_ = ExecutionEngine::createJIT(mod, &errmsg);
        if (!llvm_engine_) {
            std::cerr << "failed to create LLVM execution engine: "
                      << errmsg << "\n";
            exit(1);
        }
    } else {
        llvm_engine_->addModule(mod);
    }
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
