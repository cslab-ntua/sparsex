/*
 * Copyright (C) 2009-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2014, Vasileios Karakasis
 * Copyright (C) 2010-2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Jit.cpp
 * \brief Just In Time compilation routines
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
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
