/*
 * Copyright (C) 2017, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2017, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file JitCompiler.hpp
 * \brief Wrapper of a Clang compiler instance. Responsible for generating LLVM
 * IR code from C99 source
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_JIT_COMPILER_HPP
#define SPARSEX_INTERNALS_JIT_COMPILER_HPP

#include <iostream>
#include <cassert>
#include <boost/tokenizer.hpp>
#include <clang/Basic/Version.inc>
#include <clang/Basic/DiagnosticOptions.h>
#include <clang/Basic/LangOptions.h>
#include <clang/Basic/TargetOptions.h>
#include <clang/Basic/TargetInfo.h>
#include <clang/CodeGen/CodeGenAction.h>
#include <clang/Frontend/FrontendOptions.h>
#include <clang/Frontend/CodeGenOptions.h>
#include <clang/Frontend/CompilerInstance.h>
#include <clang/Frontend/CompilerInvocation.h>
#include <clang/Frontend/TextDiagnosticBuffer.h>
#include <clang/Frontend/TextDiagnosticPrinter.h>
#include <clang/Lex/PreprocessorOptions.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Module.h>

#include <sparsex/internals/JitConfig.hpp>
#include <sparsex/internals/JitUtil.hpp>

using namespace std;

namespace sparsex {
  namespace jit {

    class ClangCompiler
    {
    public:
      enum Options {
	IncludePathSystem,
	IncludePathUser
      };
    
      ClangCompiler();
      ~ClangCompiler() {}
    
      unique_ptr<llvm::Module> Compile(const string &source) const;

      void KeepTemporaries(bool keep)
      {
	keep_temporaries_ = keep;
      }

      void SetLogStream(ostream *log)
      {
	assert(log && "passed ostream is NULL");
      }
    
      bool DebugMode()
      {
	return debug_mode_;
      }

      void SetDebugMode(bool debug)
      {
	debug_mode_ = debug;
      }

      void AddIncludeSearchPath(const string &inc_path, Options type);

    private:
      shared_ptr<clang::CompilerInvocation> invocation_;
      unique_ptr<clang::CompilerInstance> instance_;
      bool keep_temporaries_;
      bool debug_mode_;
    };
  
  } // end of namespace jit
} // end of namespace sparsex  

#endif // SPARSEX_INTERNALS_JIT_COMPILER_HPP
