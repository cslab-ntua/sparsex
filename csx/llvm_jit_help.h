#ifndef LLVM_JIT_HELP_H__
#define LLVM_JIT_HELP_H__

#include "llvm/ADT/StringMap.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Module.h"
#include "llvm/BasicBlock.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ModuleProvider.h"


#include <map>

using namespace llvm;

IRBuilder<> *llvm_hook_builder(Module *M, const char *name);

BasicBlock *llvm_hook_newbb(Module *M, const char *name, BasicBlock **BBnext);
BasicBlock *llvm_hook_newbb(Module *M, const char *name, Function *Parent, BasicBlock **BBnext);

Module *ModuleFromFile(const char *file);
void ModuleToFile(Module *M, const char *file);

void LinkFileToModule(Module *Module, const char *file);
void doOptimize(Module *M);


Function *doCloneFunction(Module *M, Function *Fn, const char *newFnName);
Function *doCloneFunction(Module *M, const char *FnName, const char *newFnName);

Function *CloneAndReplaceHook(Module *M, Function *CallerFn, Function *HookReplaceFn,
                         const char *HookName, const char *newFnName);

bool InlineFntoFn(Function *Callee, Function *Caller);

#if 0
unsigned long getFnSize(const char *FnName);
void printSymTable();
#endif

void CheckPromoteAllocas(Value *V);

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

class SingleModule {
public:
	typedef std::map<const char *, Module *> ModMap;
	typedef std::map<Module *, unsigned int> RefMap;
	typedef std::map<Module *, ExecutionEngine *> JitMap;

	static ModMap modules;
	static RefMap refs;
	static JitMap jits;

	static Module *getM(const char *MName);
	static ExecutionEngine *getJIT(Module *M);
};

#endif /* LLVM_JIT_HELP_H__ */
