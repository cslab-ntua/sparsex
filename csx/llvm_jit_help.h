/* -*- C++ -*-
 *
 * llvm_jit_help.h -- JIT helper functions.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LLVM_JIT_HELP_H__
#define LLVM_JIT_HELP_H__

#include "llvm/ADT/StringMap.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Module.h"
#include "llvm/BasicBlock.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/Analysis/Verifier.h"

#include <map>

using namespace llvm;

/*
 * LLVM  helper functions
 */

/**
 * Utility functions:
 *  simple wrappers for common LLVM operations
 */

// load an llvm module from a bitcode file
Module *ModuleFromFile(const char *file, LLVMContext &Ctx);
// store an llvm module to a bitcode file
void ModuleToFile(Module *M, const char *file);
// link code from a bitcode file to an llvm module
void LinkFileToModule(Module *Module, const char *file, LLVMContext &Ctx);
// do optimization passes
void doOptimize(Module *M);
// helper functions for inlining
bool InlineAndRemoveFn(Function *Fn);
bool InlineFntoFn(Function *Callee, Function *Caller);
// wrappers for llvm's, CloneFunction()
Function *doCloneFunction(Module *M, const char *FnName, const char *newFnName);
Function *doCloneFunction(Module *M, Function *Fn, const char *newFnName);
Function *CloneAndReplaceFn(Module *M,
                            Function *CallerFn,
                            Function *ReplaceFn,
                            const char *FnName,
                            const char *newFnName);
// crete a JIT execution engine from a Module
ExecutionEngine *mkJIT(Module *M);

/**
 * Hook functions:
 *  Instead of building all LLVM code from the start, we use C templates,
 *  which we compile in LLVM bytecode. We use dummy functions (called
 *  hooks) to mark places in these templates.
 */

// This function replaces a hook function with two basic blocks,
// for the purpose of replacing the hook with other instructions.
//  - The returned BB is where the new instructions should be placed,
//  - BBnext is where control flow should continue, after the new instructions
//    have finished
BasicBlock *llvm_hook_newbb(Module *M, const char *name, BasicBlock **BBnext);
// same as before, but hook should be insinde Parent
BasicBlock *llvm_hook_newbb(Module *M, const char *name, Function *Parent, BasicBlock **BBnext);



/**
 * Annotations:
 *  Annotations are used in llvm template sto mark variables, so
 *  that we can detect them later in the bytecode.
 */

/**
 * Class for managing annotations (llvm.var.annotation)
 *  - find annotations in a module (annotations map strings to LLVM values)
 *  - query values, based on strings
 */
class Annotations {
    typedef StringMap<Value *> AnnotMap;
    // map from strings to LLVM values
    AnnotMap *map;

public:
    Annotations() : map(NULL) {}

    // update with annotations from M,
    // if Parent is set, consider values only inside parent
    void update(Module *M, Function *Parent=NULL);

    // get a value, for a given annotation (NULL if annotation doesnt exist)

    Value *getValue(const char *annotation);
    // dump annotations to std::cout
    void dump();
};


#endif /* LLVM_JIT_HELP_H__ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
