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
#include "template_text.h"

#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/Module.h>
#include <llvm/LLVMContext.h>


using namespace llvm;

namespace csx {

/**
 *  Singleton wrapper to an LLVM execution engine.
 */
class CsxExecutionEngine {
public:
    static CsxExecutionEngine &CreateEngine()
    {
        static CsxExecutionEngine engine;
        return engine;
    }

    void AddModule(Module *mod);

    void RemoveModule(Module *mod)
    {
        llvm_engine_->removeModule(mod);
    }

    void *GetPointerToFunction(Function *fn)
    {
        return llvm_engine_->getPointerToFunction(fn);
    }

    ~CsxExecutionEngine()
    {
        delete llvm_engine_;
    };

private:
    CsxExecutionEngine() { };
    CsxExecutionEngine(const CsxExecutionEngine &); // Do not implement
    void operator=(const CsxExecutionEngine &);     // Do not implement

    //static CsxExecutionEngine *engine_;
    ExecutionEngine *llvm_engine_;
};

/// Initialize code generation. Call it once before start generating code
CsxExecutionEngine &CsxJitInit(void);

/**
 *  Code generator for a CSX matrix (one per thread)
 */
class CsxJit {

public:
    CsxJit(CsxManager *csxmg, CsxExecutionEngine *engine, unsigned int tid = 0);

    ~CsxJit() {
        delete context_;
        delete compiler_;
    };

    void GenCode(std::ostream &log);
    void *GetSpmvFn() const;

private:
    // Compile C99 source code into an LLVM module
    Module *DoCompile(const std::string &source) const;

    //
    // Obsolete -- all code optimizations are handled by ClangCompiler
    // 
    void DoOptimizeModule();

    void DoNewRowHook(std::map<std::string, std::string> &hooks,
                      std::ostream &log, bool rowjmp) const;
    void DoSpmvFnHook(std::map<std::string, std::string> &hooks,
                      std::ostream &log);

    //
    //  Code generation functions for each pattern type
    // 
    std::string DoGenDeltaCase(int delta_bits);
    std::string DoGenLinearCase(SpmIterOrder type, int delta);
    std::string DoGenBlockCase(SpmIterOrder type, int r, int c);
    TemplateText *GetMultTemplate(SpmIterOrder type);

    CsxManager  *csxmg_;
    Module      *module_;
    CsxExecutionEngine *engine_;
    LLVMContext *context_;
    ClangCompiler *compiler_;
    unsigned int thread_id_;
    std::map<SpmIterOrder, TemplateText *> mult_templates_;
};

} // end csx namespace

#endif // CSX_JIT_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
