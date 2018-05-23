/*
 * Copyright (C) 2017, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2017, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file JitCompiler.cpp
 * \brief JIT compiler using Clang.
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2017
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/JitCompiler.hpp>

using namespace std;

namespace sparsex {
  namespace jit {

    ClangCompiler::ClangCompiler()
      : invocation_(),
	instance_(),
	keep_temporaries_(false),
	debug_mode_(false)
    {
      // Create a new CompilerInstance
      instance_.reset(new clang::CompilerInstance());
      // Create a DiagnosticsEngine
      instance_->createDiagnostics();
      if (!instance_->hasDiagnostics()) {
	cerr << "[ERROR]: createDiagnostics() failed" << endl;
	exit(1);
      }
      // Set up the diagnostic buffer for reporting errors
      instance_->getDiagnostics().setClient(new clang::TextDiagnosticBuffer());
      // Set target
      shared_ptr<clang::TargetOptions> to = make_shared<clang::TargetOptions>();
      to->Triple = llvm::sys::getDefaultTargetTriple();
      instance_->setTarget(clang::TargetInfo::CreateTargetInfo
			   (instance_->getDiagnostics(), to));
      assert(instance_->hasTarget());
      // Create a CompilerInvocation
      invocation_.reset(new clang::CompilerInvocation());
      const char *const dummy_argv[] = { "" };
      clang::CompilerInvocation::CreateFromArgs(*invocation_,
						dummy_argv,
						dummy_argv,
						instance_->getDiagnostics());
      // All options of the compiler instance should be set through the
      // invocation
#if CLANG_VERSION_MAJOR == 5
      invocation_->setLangDefaults(*invocation_->getLangOpts(),
				   clang::InputKind::C,
				   instance_->getTarget().getTriple(),
				   instance_->getPreprocessorOpts(),
				   clang::LangStandard::Kind::lang_c99);
#elif CLANG_VERSION_MAJOR == 4
      invocation_->setLangDefaults(*invocation_->getLangOpts(),
				   clang::InputKind::IK_C,
				   instance_->getTarget().getTriple(),
				   instance_->getPreprocessorOpts(),
				   clang::LangStandard::lang_c99);
#endif
      clang::LangOptions &lang_opts = *invocation_->getLangOpts();
      lang_opts.EmitAllDecls = true;
      // Setup preprocessor options
      clang::PreprocessorOptions &preproc_opts =
	invocation_->getPreprocessorOpts();
      if (debug_mode_) {
	preproc_opts.addMacroDef("SPX_DEBUG");
      } else {
	preproc_opts.addMacroDef("NDEBUG");
      }
      // Setup target options
      invocation_->getTargetOpts().Triple = llvm::sys::getDefaultTargetTriple();
      // Setup diagnostic options
      clang::DiagnosticOptions &diag_opts = invocation_->getDiagnosticOpts();
      diag_opts.Warnings.push_back("all");     // -Wall
      diag_opts.Pedantic = 1;                  // -pedantic
      diag_opts.ShowColors = 1;                // be fancy ;)
      // Setup header search options
      // Add a user-defined include path, if specified
      char *user_inc_path = getenv("SPX_JIT_INC_PATH");
      if (user_inc_path)
	AddIncludeSearchPath(user_inc_path, IncludePathSystem);
      AddIncludeSearchPath(SPX_JIT_INCLUDE, IncludePathSystem);
      AddIncludeSearchPath(CLANG_INC_SEARCH_PATH, IncludePathSystem);
      // Setup code generation options
      clang::CodeGenOptions &codegen_opts = invocation_->getCodeGenOpts();
      codegen_opts.OptimizationLevel = (debug_mode_) ? 0 : 3;
      codegen_opts.setInlining(clang::CodeGenOptions::NormalInlining);
      codegen_opts.DataSections = debug_mode_;
      codegen_opts.FunctionSections = debug_mode_;
      codegen_opts.EmitDeclMetadata = debug_mode_;
      codegen_opts.UnrollLoops = !debug_mode_;
      codegen_opts.RelaxedAliasing = !debug_mode_;
      // Assign this invocation to the instance
      instance_->setInvocation(invocation_);
      // Create FileManager
      instance_->createFileManager();
      // Create SourceManager
      instance_->createSourceManager(instance_->getFileManager());
      // Create Preprocessor using the invocation, file, and source managers
      instance_->createPreprocessor(clang::TU_Complete);
    }

    unique_ptr<llvm::Module>
    ClangCompiler::Compile(const string &source) const
    {
      // Write the source to a temporary file
      string temp_tmpl = ".tmp_XXXXXX";
      const char *tmpfile = UniqueFilename(temp_tmpl);
      SourceToFile(tmpfile, source);

      // Set the input file for compilation
      clang::FrontendOptions &frontend_opts = instance_->getFrontendOpts();
      frontend_opts.Inputs.clear(); // clear any old inputs
#if CLANG_VERSION_MAJOR == 5 || CLANG_VERSION_MAJOR == 6
      frontend_opts.Inputs.push_back(clang::FrontendInputFile
				     (tmpfile, clang::InputKind::C));
#elif CLANG_VERSION_MAJOR == 4
      frontend_opts.Inputs.push_back(clang::FrontendInputFile
				     (tmpfile, clang::InputKind::IK_C));
#endif
      // Compile and emit LLVM IR
      llvm::LLVMContext *context = new llvm::LLVMContext();
      unique_ptr<clang::CodeGenAction> llvm_codegen
	(new clang::EmitLLVMOnlyAction(context));
      if (!instance_->ExecuteAction(*llvm_codegen)) {
	cerr << "[ERROR]: compilation failed, "
	     << "generated source is in " << tmpfile << "\n";
	if (!keep_temporaries_)
	  RemoveFile(tmpfile);
	exit(1);
      }

      // Remove input file
      if (!keep_temporaries_ && tmpfile)
	RemoveFile(tmpfile);

      // Return the compiled module
      // Note: takeModule() returns a unique_ptr<>
      return llvm_codegen->takeModule();
    }

    void
    ClangCompiler::AddIncludeSearchPath(const string &inc_path, Options type)
    {
      clang::frontend::IncludeDirGroup inc_group = clang::frontend::Angled;
      if (type == IncludePathSystem) {
	inc_group = clang::frontend::System;
      }

      clang::HeaderSearchOptions &header_search =
	invocation_->getHeaderSearchOpts();

      boost::char_separator<char> path_sep(":");
      boost::tokenizer<boost::char_separator<char> > path_tokens(inc_path,
								 path_sep);
      boost::tokenizer<boost::char_separator<char> >::iterator tok_iter;
      for (tok_iter = path_tokens.begin();
	   tok_iter != path_tokens.end(); ++tok_iter)
	header_search.AddPath(*tok_iter, inc_group, false, false);
    }

  } // end of namespace jit
} // end of namespace sparsex
