#ifndef LLVM_JIT_HELP_H__
#define LLVM_JIT_HELP_H__

#include "llvm/ADT/StringMap.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Module.h"
#include "llvm/BasicBlock.h"

using namespace llvm;

IRBuilder<> *llvm_hook_builder(Module *M, const char *name);
BasicBlock *llvm_hook_newbb(Module *M, const char *name, BasicBlock **BBnext);

Module *ModuleFromFile(const char *file);
void LinkFileToModule(Module *Module, const char *file);
void doOptimize(Module *M);
Function *CloneAndReplaceHook(Module *M, Function *CallerFn, Function *HookReplaceFn,
                         const char *HookName, const char *newFnName);

bool InlineFntoFn(Function *Callee, Function *Caller);

#if 0
unsigned long getFnSize(const char *FnName);
void printSymTable();
#endif

typedef StringMap<Value *> AnnotationMap;

// Class for accessing annotations (llvm.var.annotation)
class Annotations {
	AnnotationMap *map;

public:
	Annotations() : map(NULL) {}

	void update(Module *M);
	void dump();
	Value *getValue(const char *annotation);
};

#endif /* LLVM_JIT_HELP_H__ */
