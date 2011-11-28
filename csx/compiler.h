/* -*- C++ -*-
 *
 * compiler.h -- Wrapper of a Clang compiler instance. Responsible for
 *               generating LLVM IR code from C99 source.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef COMPILER_H__
#define COMPILER_H__

#include <clang/Frontend/DiagnosticOptions.h>
#include <clang/Frontend/CompilerInstance.h>
#include <llvm/Module.h>
#include <iostream>
#include <cassert>

using namespace clang;
using namespace llvm;

//
// Clang compiler installation prefix; necessary to correctly initialize the
// compiler's resources path.
// 
#ifndef CLANG_PREFIX
#   define CLANG_PREFIX "/usr/local"
#endif

class ClangCompiler
{
public:
    ClangCompiler(const char *prefix = CLANG_PREFIX);

    ~ClangCompiler() {
        // Take back the ownership of the invocation so as to avoid double
        // delete() corruption; the invocation will be released from our dtor
        compiler_->takeInvocation();
    };

    Module *Compile(const std::string &source, LLVMContext *context) const;

    void KeepTemporaries(bool keep)
    {
        keep_temporaries_ = keep;
    }

    void SetLogStream(std::ostream *log)
    {
        assert(log && "passed ostream is NULL");
        log_stream_ = log;
    }

private:
    OwningPtr<CompilerInvocation> invocation_;
    OwningPtr<CompilerInstance> compiler_;
    bool keep_temporaries_;
    std::ostream *log_stream_;
};

#endif // COMPILER_H__
