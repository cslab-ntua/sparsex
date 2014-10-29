/*
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Compiler.cpp
 * \brief Wrapper of a Clang compiler instance. Responsible for generating LLVM
 * IR code from C99 source
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/Compiler.hpp>
#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/JitUtil.hpp>
#include <clang/Basic/Diagnostic.h>
#include <clang/Basic/Version.h>
#include <clang/CodeGen/CodeGenAction.h>
#include <clang/Frontend/CodeGenOptions.h>
#include <clang/Frontend/HeaderSearchOptions.h>
#include <clang/Frontend/TextDiagnosticPrinter.h>
#include <llvm/ADT/StringRef.h>
#include <llvm/Target/TargetOptions.h>
#include <boost/tokenizer.hpp>
#include <fstream>

using namespace boost;
using namespace std;

namespace sparsex {
namespace jit {

ClangCompiler::ClangCompiler()
    : invocation_(new CompilerInvocation()),
      compiler_(new CompilerInstance()), keep_temporaries_(false),
      debug_mode_(false), log_stream_(&cerr)
{
    // Set-up the clang compiler
    TextDiagnosticPrinter *diag_client =
        new TextDiagnosticPrinter(errs(), DiagnosticOptions());

    IntrusiveRefCntPtr<DiagnosticIDs> diag_id(new DiagnosticIDs());
    DiagnosticsEngine diags(diag_id, diag_client);

    // Create a dummy invocation of the compiler, so as to set it up
    // and we will replace the source file afterwards
    const char *const dummy_argv[] = { "" };
    CompilerInvocation::CreateFromArgs(*invocation_, dummy_argv, dummy_argv,
                                       diags);
    // Compile C99
    invocation_->setLangDefaults(IK_C, LangStandard::lang_c99);

    // Setup the include path
    AddIncludeSearchPath(CLANG_INC_SEARCH_PATH, IncludePathSystem);

    // Setup diagnostic options
    DiagnosticOptions &diag_options = invocation_->getDiagnosticOpts();

    // FIXME: The following option crashes Clang 3.0 on certain systems!
    // diag_options.Warnings.push_back("all");     // -Wall
    diag_options.Pedantic = 1;                  // -pedantic
    diag_options.ShowColors = 1;                // be fancy ;)

    // Setup code generation options
    SetCodeGenOptions();
}

Module *ClangCompiler::Compile(const string &source,
                               LLVMContext *context) const
{
    // write the source to a temporary file and invoke the compiler
    string temp_tmpl = ".tmp_XXXXXX";
    const char *tmpfile = UniqueFilename(temp_tmpl);
    //const char *tmpfile = "temp.c";

    SourceToFile(tmpfile, source);

    FrontendOptions &opts = invocation_->getFrontendOpts();
    opts.Inputs.clear();    // clear any old inputs
    opts.Inputs.push_back(make_pair(IK_C, tmpfile));

//    compiler_->setInvocation(invocation_.get());
    compiler_->setInvocation(invocation_);

    // Setup diagnostics for the compilation process itself
    const char *const dummy_argv[] = { "" };
    compiler_->createDiagnostics(1, dummy_argv);
    if (!compiler_->hasDiagnostics()) {
        *log_stream_ << "createDiagnostics() failed\n";
        return 0;
    }

    // Compile and emit LLVM IR
    OwningPtr<CodeGenAction> llvm_codegen(new EmitLLVMOnlyAction(context));
    if (!compiler_->ExecuteAction(*llvm_codegen)) {
        *log_stream_ << "compilation failed: "
                     << "generated source is in " << tmpfile << "\n";
        return 0;
    }

    // Remove input file and return the compiled module
    if (!keep_temporaries_)
        RemoveFile(tmpfile);
    return llvm_codegen->takeModule();
}

void ClangCompiler::AddIncludeSearchPath(const string &inc_path, Options type)
{
    bool is_user_path = true;
    frontend::IncludeDirGroup inc_group = frontend::Angled;
    if (type == IncludePathSystem) {
        is_user_path = false;
        inc_group = frontend::System;
    }

    HeaderSearchOptions &header_search =
        invocation_->getHeaderSearchOpts();

    char_separator<char> path_sep(":");
    tokenizer<char_separator<char> > path_tokens(inc_path, path_sep);
    tokenizer<char_separator<char> >::iterator tok_iter;
    for (tok_iter = path_tokens.begin();
         tok_iter != path_tokens.end(); ++tok_iter)
        header_search.AddPath(*tok_iter, inc_group,
                              is_user_path, false, false);
}

void ClangCompiler::SetCodeGenOptions()
{
    CodeGenOptions &codegen_options = invocation_->getCodeGenOpts();
    codegen_options.OptimizationLevel = (debug_mode_) ? 0 : 3;
    codegen_options.DebugInfo = debug_mode_;
    codegen_options.DataSections = debug_mode_;
    codegen_options.FunctionSections = debug_mode_;
    codegen_options.EmitDeclMetadata = debug_mode_;
    codegen_options.UnrollLoops = !debug_mode_;
    codegen_options.RelaxedAliasing = !debug_mode_;
    codegen_options.Inlining = (debug_mode_) ? CodeGenOptions::NoInlining :
        CodeGenOptions::NormalInlining;
    codegen_options.RelaxedAliasing = !debug_mode_;
    // Enable debugging from GDB
    JITEmitDebugInfo = debug_mode_; // llvm global

    PreprocessorOptions &preproc_options = invocation_->getPreprocessorOpts();
    if (debug_mode_) {
        preproc_options.addMacroDef("SPX_DEBUG");
    } else {
        preproc_options.addMacroDef("NDEBUG");
    }
}

} // end of namespace jit
} // end of namespace sparsex
