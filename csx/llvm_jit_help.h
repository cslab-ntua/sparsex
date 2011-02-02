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

// Some LLVM  helper functions

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

ExecutionEngine *mkJIT(Module *M);


/**
 * Annotations:
 */
typedef StringMap<Value *> AnnotationMap;

// Class for accessing annotations (llvm.var.annotation)
class Annotations {
    AnnotationMap *map;

public:
    Annotations() : map(NULL) {}

    void update(Module *M, Function *Parent=NULL);
    void dump();
    Value *getValue(const char *annotation);
};


#endif /* LLVM_JIT_HELP_H__ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
