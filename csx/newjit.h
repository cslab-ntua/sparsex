/* -*- C++ -*-
 *
 * jit.h -- Just In Time compilation routines.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2011, Vasileios Karakasis
 * Copyright (C) 2010-2011, Theodors Goudouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_JIT_H__
#define CSX_JIT_H__

#include "csx.h"
#include "compiler.h"

#include <llvm/Module.h>
#include <llvm/LLVMContext.h>

/// Initialize code generation. Call it once before start generating code
void CsxJitInit(void);

using namespace llvm;

namespace csx {

/**
 *  Code generator for a CSX matrix (one-per-thread)
 */
class CsxJit {

public:
    CsxJit(CsxManager *csxmg, unsigned int id = 0);

    ~CsxJit() {
        delete module_;
        delete context_;
        delete compiler_;
    };

    void GenCode(std::ostream &log);
    void *GetSpmvFn() const;

private:
    // Compile C99 source code into an LLVM module
    Module *DoCompile(const std::string &source) const;

    void DoOptimize(Module *module) const;

    void DoNewRowHook(std::map<std::string, std::string> &hooks,
                      std::ostream &log) const;
    void DoSpmvFnHook(std::map<std::string, std::string> &hooks,
                      std::ostream &log) const;

    //
    //  Code generation functions for each pattern type
    // 
    void DoGenDeltaCase(int delta_bytes, std::string &code) const;
    void DoGenHorizCase(int delta, std::string &code) const;
    void DoGenVertCase(int delta, std::string &code) const;
    void DoGenDiagCase(int delta, std::string &code) const;
    void DoGenRDiagCase(int delta, std::string &code) const;
    void DoGenBlockRowCase(int r, int c, std::string &code) const;
    void DoGenBlockColCase(int r, int c, std::string &code) const;

    CsxManager  *csxmg_;
    Module      *module_;
    LLVMContext *context_;
    ClangCompiler *compiler_;
};

} // end csx namespace

#endif // CSX_JIT_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
