#include <iostream>
#include <fstream>
#include <cstring>

#include "llvm/Module.h"
#include "llvm/Type.h"
#include "llvm/Function.h"
#include "llvm/DerivedTypes.h"
#include "llvm/ValueSymbolTable.h"
#include "llvm/TypeSymbolTable.h"
#include "llvm/Constant.h"
#include "llvm/Instruction.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/BasicBlock.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Bitcode/ReaderWriter.h"
#include "llvm/Linker.h"
#include "llvm/Target/TargetData.h"
#include "llvm/GlobalValue.h" /* Linkage */
#include "llvm/PassManager.h"
#include "llvm/Transforms/IPO.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/Verifier.h"
#include "llvm/Transforms/Utils/Cloning.h" /* InlineFunction */


#include "llvm_jit_help.h"

using namespace llvm;

Module *ModuleFromFile(const char *file)
{
	std::string Error;
	MemoryBuffer *MB = MemoryBuffer::getFile(file, &Error);
	if (!MB){
		std::cerr << " MemoryBuffer::getFile " << Error << " " << file << "\n";
		exit(1);
	}
	Module *M = ParseBitcodeFile(MB, &Error);
	if (!M){
		std::cerr << "ParseBitCodeFile:" << Error << "\n";
		exit(1);
	}

	return M;
}

void ModuleToFile(Module *M, const char *file)
{
	std::string Error;

	std::ofstream out;
	out.open(file);
	WriteBitcodeToFile(M, out);
	out.close();
}

void LinkFileToModule(Module *M, const char *file)
{
	std::string Err;
	Module *newM = ModuleFromFile(file);
	int err = Linker::LinkModules(M, newM, &Err);
	if (err){
		std::cerr << "Error in linking " << file << " " << Err << "\n";
		exit(1);
	}
	delete newM;
}

/* check that the function is void fn(void) */
static inline int is_type_ok(const Type *Ty)
{
	const FunctionType *FnTy = cast<FunctionType>(Ty);
	return (FnTy->getNumParams() == 0) &&
	       (FnTy->getReturnType() == Type::VoidTy);
}


static Instruction *get_hook(Module *M, const char *name)
{
	Function *Fn = M->getFunction(name);

	const Type *FnTy = Fn->getFunctionType();
	if (!is_type_ok(FnTy)){
		std::cout << "Error in Function Type\n";
		return NULL;
	}

	/* get the hook, and verify it is used only once */
	Value::use_iterator I = Fn->use_begin(), E = Fn->use_end();
	if (I == E){
		std::cout << "No uses of the hook\n";
		return NULL;
	}
	Instruction *Hook = cast<Instruction>(*I);
	if (++I != E) {
		std::cout << "More than one uses of the hook\n";
		return NULL;
	}

	return Hook;
}

// Get the first call instruction of Fn in CallerFn
CallInst *getCallInFn(Function *Fn, Function *CallerFn)
{
	CallInst *CI;
	Value::use_iterator I = Fn->use_begin(), E = Fn->use_end();
	for ( ; I != E; I++){
		CI = cast<CallInst>(*I);
		if (CI->getParent()->getParent() == CallerFn){
			return CI;
		}
	}

	std::cout << "Function not found\n";
	return NULL;
}


Function *CloneAndReplaceHook(Module *M, Function *CallerFn, Function *HookReplaceFn,
                         const char *HookName, const char *newFnName)
{
	// get Hook Function
	Function *HookFn = M->getFunction(HookName);

	assert(HookFn);
	//printf("Uses of %s: %u\n", HookFn->getNameStart(), HookFn->getNumUses());
	// Clone Caller Function and put it into module
	Function *newFn = CloneFunction(CallerFn);
	newFn->setName(newFnName);
	M->getFunctionList().push_back(newFn);

	// get call instruction to replace in new function
	//printf("Uses of %s: %u\n", HookFn->getNameStart(), HookFn->getNumUses());
	CallInst *CallInstr = getCallInFn(HookFn, newFn);
	CallInstr->setOperand(0, HookReplaceFn);
	//printf("Uses of %s: %u\n", HookFn->getNameStart(), HookFn->getNumUses());
	verifyModule(*M, AbortProcessAction, 0);

	return newFn;
}

#if 0
IRBuilder *llvm_hook_builder(Module *M, const char *name)
{
	Instruction *Hook = get_hook(M,name);

	/* remove the hook and create an IRBuilder in the proper place */
	BasicBlock *BB = Hook->getParent();
	BasicBlock::iterator BI(Hook);
	BI = BB->getInstList().erase(BI);
	IRBuilder *ret = new IRBuilder(BB, BI);
	return ret;
}
#endif

// This function replaces a hook function with two basic blocks,
// for the purpose of replacing the hook with other instructions.
// The return basic  block is where the new instructions should be placed,
// while BBnext, is where control flow should continue
BasicBlock *llvm_hook_newbb(Module *M, const char *name, BasicBlock **BBnext)
{
	Instruction *Hook = get_hook(M,name);

	BasicBlock *BB = Hook->getParent();
	BasicBlock::iterator BI(Hook);
	BI = BB->getInstList().erase(BI);

	BasicBlock *BB0 = BB->splitBasicBlock(BI, "split_start");
	BasicBlock *BB1 = BB0->splitBasicBlock(BB0->begin(), "split_end");

	/* BB0 should have only one (uncond) branch instruction. remove it */
	BB0->getInstList().erase(BB0->begin());
	assert(BB0->begin() == BB0->end());

	*BBnext = BB1;
	return BB0;
}

/* XXX: Copied from JITEmitter.cpp  */

/// JitSymbolEntry - Each function that is JIT compiled results in one of these
/// being added to an array of symbols.  This indicates the name of the function
/// as well as the address range it occupies.  This allows the client to map
/// from a PC value to the name of the function.
struct JitSymbolEntry {
  const char *FnName;   // FnName - a strdup'd string.
  void *FnStart;
  intptr_t FnSize;
};


struct JitSymbolTable {
  /// NextPtr - This forms a linked list of JitSymbolTable entries.  This
  /// pointer is not used right now, but might be used in the future.  Consider
  /// it reserved for future use.
  JitSymbolTable *NextPtr;

  /// Symbols - This is an array of JitSymbolEntry entries.  Only the first
  /// 'NumSymbols' symbols are valid.
  JitSymbolEntry *Symbols;

  /// NumSymbols - This indicates the number entries in the Symbols array that
  /// are valid.
  unsigned NumSymbols;

  /// NumAllocated - This indicates the amount of space we have in the Symbols
  /// array.  This is a private field that should not be read by external tools.
  unsigned NumAllocated;
};

#if 0
extern JitSymbolTable *__jitSymbolTable;
#endif

unsigned long getFnSize(const char *FnName)
{
	#if 0
	unsigned int i;
	unsigned long len0 = strlen(FnName);
	for (i=0; i<__jitSymbolTable->NumSymbols; i++){
		JitSymbolEntry *entry = &__jitSymbolTable->Symbols[i];
		unsigned long len1 = strlen(entry->FnName);
		if (len0 != len1)
			continue;
		if ( memcmp(FnName, entry->FnName, len0) == 0 ){
			return entry->FnSize;
		}
	}
	#endif

	return 0;
}

void printSymTable()
{
	#if 0
	unsigned int i;
	for (i=0; i<__jitSymbolTable->NumSymbols; i++){
		JitSymbolEntry *entry = &__jitSymbolTable->Symbols[i];
		std::cout << entry->FnName << " : " << entry->FnStart << " " << entry->FnSize << "\n";
	}
	#endif
}


inline void addPass(PassManager &PM, Pass *P) {
	PM.add(P);
}

void doGlobalDCE(Module *M)
{
	PassManager PM;
	PM.add(new TargetData(M));
	addPass(PM, createGlobalDCEPass());
	PM.run(*M);
}

void doOptimize(Module *M)
{
	/* copied from opt.cpp */

	PassManager PM;
	PM.add(new TargetData(M));

	addPass(PM, createLowerSetJmpPass());          // Lower llvm.setjmp/.longjmp

	// If the -strip-debug command line option was specified, do it.
	//if (StripDebug)
	//addPass(PM, createStripSymbolsPass(true));
	//if (DisableOptimizations) return;

	addPass(PM, createRaiseAllocationsPass());     // call %malloc -> malloc inst
	addPass(PM, createCFGSimplificationPass());    // Clean up disgusting code
	addPass(PM, createPromoteMemoryToRegisterPass());// Kill useless allocas
	addPass(PM, createGlobalOptimizerPass());      // Optimize out global vars
	addPass(PM, createGlobalDCEPass());            // Remove unused fns and globs
	addPass(PM, createIPConstantPropagationPass());// IP Constant Propagation
	addPass(PM, createDeadArgEliminationPass());   // Dead argument elimination
	addPass(PM, createInstructionCombiningPass()); // Clean up after IPCP & DAE
	addPass(PM, createCFGSimplificationPass());    // Clean up after IPCP & DAE
	addPass(PM, createPruneEHPass());              // Remove dead EH info

	//if (!DisableInline)
	addPass(PM, createFunctionInliningPass());   // Inline small functions
	addPass(PM, createArgumentPromotionPass());    // Scalarize uninlined fn args

	addPass(PM, createTailDuplicationPass());      // Simplify cfg by copying code
	addPass(PM, createSimplifyLibCallsPass());     // Library Call Optimizations
	addPass(PM, createInstructionCombiningPass()); // Cleanup for scalarrepl.
	addPass(PM, createJumpThreadingPass());        // Thread jumps.
	addPass(PM, createCFGSimplificationPass());    // Merge & remove BBs
	addPass(PM, createScalarReplAggregatesPass()); // Break up aggregate allocas
	addPass(PM, createInstructionCombiningPass()); // Combine silly seq's
	addPass(PM, createCondPropagationPass());      // Propagate conditionals

	addPass(PM, createTailCallEliminationPass());  // Eliminate tail calls
	addPass(PM, createCFGSimplificationPass());    // Merge & remove BBs
	addPass(PM, createReassociatePass());          // Reassociate expressions
	addPass(PM, createLoopRotatePass());
	addPass(PM, createLICMPass());                 // Hoist loop invariants
	addPass(PM, createLoopUnswitchPass());         // Unswitch loops.
	addPass(PM, createLoopIndexSplitPass());       // Index split loops.
	addPass(PM, createInstructionCombiningPass()); // Clean up after LICM/reassoc
	addPass(PM, createIndVarSimplifyPass());       // Canonicalize indvars
	addPass(PM, createLoopUnrollPass());           // Unroll small loops
	addPass(PM, createInstructionCombiningPass()); // Clean up after the unroller
	addPass(PM, createGVNPass());                  // Remove redundancies
	addPass(PM, createMemCpyOptPass());            // Remove memcpy / form memset
	addPass(PM, createSCCPPass());                 // Constant prop with SCCP

	// Run instcombine after redundancy elimination to exploit opportunities
	// opened up by them.
	addPass(PM, createInstructionCombiningPass());
	addPass(PM, createCondPropagationPass());      // Propagate conditionals

	addPass(PM, createDeadStoreEliminationPass()); // Delete dead stores
	addPass(PM, createAggressiveDCEPass());        // SSA based 'Aggressive DCE'
	addPass(PM, createCFGSimplificationPass());    // Merge & remove BBs
	addPass(PM, createStripDeadPrototypesPass());  // Get rid of dead prototypes
	addPass(PM, createDeadTypeEliminationPass());  // Eliminate dead types
	addPass(PM, createConstantMergePass());        // Merge dup global constants

	PM.run(*M);
}

// Force inline and try to remove a function by setting intenal linkage
//   - Set the linkage of a function to internal
//   - Inline its uses
//   - call a dead code elimination pass to remove its uses (hopefully)
//   - return true if removal succeeded
bool InlineAndRemoveFn(Function *Fn)
{
	Fn->setLinkage(GlobalValue::InternalLinkage);
	Value::use_iterator ui = Fn->use_begin();
	Value::use_iterator ue = Fn->use_end();
	for (; ui != ue; ){
		CallInst *CI = cast<CallInst>(*ui++);
		InlineFunction(CI);
	}
	std::string FnName = Fn->getName();
	Module *M = Fn->getParent();
	doGlobalDCE(M);

	return (M->getFunction(FnName) != NULL);
}

// Force inline of a callee  to the first call instance of caller
// return false if call wasn't found
bool InlineFntoFn(Function *Callee, Function *Caller)
{
	Value::use_iterator ui = Callee->use_begin();
	Value::use_iterator ue = Callee->use_end();
	for (; ui != ue; ){
		CallInst *CI = cast<CallInst>(*ui++);
		if (CI->getParent()->getParent() == Caller){
			if ( !InlineFunction(CI) ){
				std::cout << __FUNCTION__ << ": Function not inlined\n";
				return false;
			}
			return true;
		}
	}

	std::cout << __FUNCTION__ << ": Function not found\n";
	return false;
}

// update map, with annotations from Module M
void Annotations::update(Module *M)
{
	Function *Fn;
	Value::use_iterator ui;

	Fn = M->getFunction("llvm.var.annotation");
	if (!Fn){
		std::cerr << "no annotations found" << std::endl;
		return;
	}

	if (this->map == NULL)
		this->map = new AnnotationMap();

	for (ui = Fn->use_begin(); ui != Fn->use_end(); ui++){
		CallInst *CI;
		Value *val;
		GlobalVariable *ptr;
		std::string str;
		// jump trough the hoops
		CI = cast<CallInst>(*ui);
		val = CI->getOperand(1)->stripPointerCasts();
		ptr = cast<GlobalVariable>(CI->getOperand(2)->stripPointerCasts());
		str = cast<ConstantArray>(ptr->getInitializer())->getAsString();
		// insert (the C string)
		(*this->map)[str.c_str()] = val;
	}
}

void Annotations::dump()
{
	AnnotationMap::iterator am = this->map->begin();

	for (; am != this->map->end(); ++am){
		std::cout << am->getKeyData() << "\n";
	}
}

Value *Annotations::getValue(const char *annotation)
{
	AnnotationMap::iterator am = this->map->find(annotation);
	assert(am != this->map->end());
	return am->getValue();
}
