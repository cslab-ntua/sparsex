/*
 * \file Jit.cpp
 *
 * \brief Just In Time compilation utilities
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * Copyright (C) 2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/Jit.hpp>

using namespace std;
using namespace llvm;

namespace sparsex {
namespace jit {

CsxExecutionEngine &CsxJitInit(void)
{
    InitializeNativeTarget();
    return CsxExecutionEngine::CreateEngine();
}

void CsxExecutionEngine::AddModule(Module *mod)
{
    if (!llvm_engine_) {
        string errmsg;
        llvm_engine_ = ExecutionEngine::createJIT(mod, &errmsg);
        if (!llvm_engine_) {
            cerr << "failed to create LLVM execution engine: "
                 << errmsg << "\n";
            exit(1);
        }
    } else {
        llvm_engine_->addModule(mod);
    }
}

} // end of namespace jit
} // end of namespace sparsex

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
