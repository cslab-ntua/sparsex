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

#ifndef SPARSEX_INTERNALS_CODE_EXECUTOR_HPP
#define SPARSEX_INTERNALS_CODE_EXECUTOR_HPP

#include <iostream>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/ExecutionEngine/MCJIT.h>
#include <llvm/IR/Module.h>

using namespace std;

namespace sparsex {
  namespace jit {

  // Singleton wrapper to an LLVM execution engine
  class CodeExecutor {
  public:

    static CodeExecutor &GetInstance()
    {
      // C++11 mandates that the initializer for a local static variable is
      // only run once, even in the presence of concurrency. So, this code
      // is thread-safe.
      static CodeExecutor instance;
      return instance;
    }
    
    ~CodeExecutor();
    // Assumes the given module has been verified
    void AddModule(unique_ptr<llvm::Module> module);
    void RemoveModule(llvm::Module *module);
    uint64_t GetFnAddress(const string fn_name);

  private:
    llvm::EngineBuilder *builder_;
    llvm::ExecutionEngine *engine_;

    CodeExecutor();    
    CodeExecutor(const CodeExecutor &) = delete; // Do not implement
    void operator=(const CodeExecutor &) = delete; // Do not implement
  };

} // end of namespace jit
} // end of namespace sparsex  

#endif  // SPARSEX_INTERNALS_CODE_EXECUTOR_HPP
