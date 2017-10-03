/*
 * Copyright (C) 2017, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2017, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CodeExecutor.cpp
 * \brief Component that dynamically loads a compiled module in memory for
          execution based on LLVM.
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2017
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/CodeExecutor.hpp>

using namespace std;

namespace sparsex {
  namespace jit {

    CodeExecutor::CodeExecutor()
      : engine_(nullptr)
    {
      LLVMInitializeNativeTarget();
      LLVMInitializeNativeAsmPrinter();
      LLVMInitializeNativeAsmParser();
      LLVMLinkInMCJIT();
    }

    CodeExecutor::~CodeExecutor()
    {
      // check if engine contains modules delete them/remove them?
      delete engine_;
      delete builder_;
      llvm::llvm_shutdown();
    }

    void CodeExecutor::AddModule(unique_ptr<llvm::Module> module)
    {
      // The ExecutionEngine takes ownership of the module
      if (!engine_) {
	string errmsg;
	builder_ = new llvm::EngineBuilder(move(module));
	builder_->setErrorStr(&errmsg);
	builder_->setEngineKind(llvm::EngineKind::JIT);
	builder_->setVerifyModules(true);
	engine_ = builder_->create();
	if (!engine_) {
	  cerr << "[ERROR]: failed to create LLVM execution engine, "
	       << errmsg
	       << endl;
	  exit(1);
	}
      } else {
	engine_->addModule(move(module));
      }

      assert(module.get() == nullptr); 
      // Code generation is deferred here
      engine_->finalizeObject();
    }
    
    void CodeExecutor::RemoveModule(llvm::Module *module)
    {
      engine_->removeModule(module);
    }

    uint64_t CodeExecutor::GetFnAddress(const string func_name)
    {
      // Resolves name mangling for C++ code
      return engine_->getFunctionAddress(func_name);
    }
  
  } // end of namespace jit
} // end of namespace sparsex  
