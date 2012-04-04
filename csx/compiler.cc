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

#include "compiler.h"
#include "jit_util.h"

#include <clang/Basic/Version.h>
#include <clang/CodeGen/CodeGenAction.h>
#include <clang/Frontend/CodeGenOptions.h>
#include <clang/Frontend/HeaderSearchOptions.h>
#include <clang/Frontend/TextDiagnosticPrinter.h>
#include <llvm/ADT/StringRef.h>
#include "llvm/Target/TargetOptions.h"
#include <fstream>

static std::string GetClangResourcePath(const char *prefix)
{
    return std::string(prefix) + "/lib/clang/" + CLANG_VERSION_STRING +
        "/include";
}

ClangCompiler::ClangCompiler(const char *prefix)
    : invocation_(new CompilerInvocation()),
      compiler_(new CompilerInstance()), keep_temporaries_(false),
      debug_mode_(false), log_stream_(&std::cerr)
{
    // Set-up the clang compiler
    TextDiagnosticPrinter *diag_client =
        new TextDiagnosticPrinter(errs(), DiagnosticOptions());

    IntrusiveRefCntPtr<DiagnosticIDs> diag_id(new DiagnosticIDs());
    Diagnostic diags(diag_id, diag_client);

    // Create a dummy invocation of the compiler, so as to set it up
    // and we will replace the source file afterwards
    const char *const dummy_argv[] = { "" };
    CompilerInvocation::CreateFromArgs(*invocation_, dummy_argv, dummy_argv,
                                       diags);
    // Compile C99
    invocation_->setLangDefaults(IK_C, LangStandard::lang_c99);

    // Setup the include path
    HeaderSearchOptions &header_search =
        invocation_->getHeaderSearchOpts();
    header_search.AddPath(GetClangResourcePath(prefix),
                          frontend::System, false, false, false);
    // FIXME: do sth more generic
    header_search.AddPath("../lib/spm", frontend::Quoted,
                          true /* user supplied */, false, false);
    header_search.AddPath("../lib/dynarray", frontend::Angled,
                          true /* user supplied */, false, false);
    header_search.AddPath("/usr/include/x86_64-linux-gnu", frontend::System,
                          false, false, false);

    // Setup diagnostic options
    DiagnosticOptions &diag_options = invocation_->getDiagnosticOpts();
    diag_options.Warnings.push_back("all");     // -Wall
    diag_options.Pedantic = 1;                  // -pedantic
    diag_options.ShowColors = 1;                // be fancy ;)

    // Setup code generation options
    SetCodeGenOptions();
}

Module *ClangCompiler::Compile(const std::string &source,
                               LLVMContext *context) const
{
    // write the source to a temporary file and invoke the compiler
    std::string temp_tmpl = ".tmp_XXXXXX";
    const char *tmpfile = UniqueFilename(temp_tmpl);
    //const char *tmpfile = "temp.c";

    SourceToFile(tmpfile, source);

    FrontendOptions &opts = invocation_->getFrontendOpts();
    opts.Inputs.clear();    // clear any old inputs
    opts.Inputs.push_back(std::make_pair(IK_C, tmpfile));

    compiler_->setInvocation(invocation_.get());

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
        preproc_options.addMacroDef("CSX_DEBUG");
    } else {
        preproc_options.addMacroDef("NDEBUG");
    }
}
