/*
 * Compiler.hpp -- Wrapper of a Clang compiler instance. Responsible for
 *                 generating LLVM IR code from C99 source.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_COMPILER_HPP
#define SPARSEX_INTERNALS_COMPILER_HPP

#include <clang/Frontend/DiagnosticOptions.h>
#include <clang/Frontend/CompilerInstance.h>
#include <llvm/Module.h>
#include <iostream>
#include <cassert>
#include <llvm/Support/ManagedStatic.h>

using namespace clang;
using namespace llvm;
using namespace std;

class ClangCompiler
{
public:
    ClangCompiler();

    ~ClangCompiler() {
        // Take back the ownership of the invocation so as to avoid double
        // delete() corruption; the invocation will be released from our dtor
        //compiler_->takeInvocation();

        /* Just let CompilerInvocation object be released through CompilerInstance */
    };

    Module *Compile(const string &source, LLVMContext *context) const;

    void KeepTemporaries(bool keep)
    {
        keep_temporaries_ = keep;
    }

    void SetLogStream(ostream *log)
    {
        assert(log && "passed ostream is NULL");
        log_stream_ = log;
    }

    bool DebugMode()
    {
        return debug_mode_;
    }

    void SetDebugMode(bool debug)
    {
        debug_mode_ = debug;
        SetCodeGenOptions();
    }

    void AddHeaderSearchPath(const char *path)
    {
        HeaderSearchOptions &header_search =
            invocation_->getHeaderSearchOpts();

        // Add path both as a quoted and angled include
        header_search.AddPath(path, frontend::Quoted,
                              true /* user supplied */, false, false);
        header_search.AddPath(path, frontend::Angled,
                              true /* user supplied */, false, false);
    }

private:
    // Set up the code generation options depending on debug mode
    void SetCodeGenOptions();
    void SetHeaderSearchOptions(const string &);
    CompilerInvocation *invocation_;
    OwningPtr<CompilerInstance> compiler_;
    bool keep_temporaries_;
    bool debug_mode_;
    ostream *log_stream_;
};

#endif // SPARSEX_INTERNALS_COMPILER_HPP
