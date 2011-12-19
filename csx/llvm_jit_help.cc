/*
 * llvm_jit_help.cc -- JIT helper functions.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C)      2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include "llvm_jit_help.h"

#include "llvm/Module.h"
#include "llvm/Type.h"
#include "llvm/Function.h"
#include "llvm/Constant.h"
#include "llvm/Instruction.h"
#include "llvm/LLVMContext.h"
#include "llvm/BasicBlock.h"
#include "llvm/Linker.h"
#include "llvm/GlobalValue.h" /* Linkage */
#include "llvm/PassManager.h"
#include "llvm/Bitcode/ReaderWriter.h"
#include "llvm/Target/TargetData.h"
#include "llvm/Transforms/Utils/Cloning.h" /* InlineFunction */
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/Verifier.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/StandardPasses.h"
#include "llvm/Support/system_error.h"

#include <iostream>
#include <fstream>
#include <cstring>

/*
 * LLVM (v2.9) helpers for JIT compilation
 */

using namespace llvm;

Module *ModuleFromFile(const char *file, LLVMContext &Ctx)
{
    OwningPtr<MemoryBuffer> membuff;
    error_code err = MemoryBuffer::getFile(StringRef(file), membuff);

    if (err) {
        std::cerr << __FUNCTION__ << "(" __FILE__ << ":" << __LINE__ << "):"
                  << err.message() << std::endl;
        exit(1);
    }

    std::string err_msg;
    Module *mod = ParseBitcodeFile(membuff.get(), Ctx, &err_msg);
    if (!mod){
        std::cerr << __FUNCTION__ << "(" __FILE__ << ":" << __LINE__ << "):"
                  << err_msg << std::endl;
        exit(1);
    }

    return mod;
}

void ModuleToFile(Module *M, const char *file)
{
    std::string error;
    raw_fd_ostream output(file, error, raw_fd_ostream::F_Binary);

    if (!error.empty()){
        std::cerr << error << ": " << file << "\n";
        exit(1);
    }

    WriteBitcodeToFile(M, output);
}

void LinkFileToModule(Module *M, const char *file, LLVMContext &Ctx)
{
    std::string Err;
    Module *newM = ModuleFromFile(file, Ctx);
    int err = Linker::LinkModules(M, newM, &Err);
    if (err){
        std::cerr << "Error in linking " << file << " " << Err << "\n";
        exit(1);
    }
    delete newM;
}

inline void addPass(PassManager &PM, Pass *P)
{
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
    /* XXX: adapted from opt.cpp, not sure how correct or non-redundant is */

    PassManager PM;
    FunctionPassManager *FPasses;

    FPasses = new FunctionPassManager(M);
    FPasses->add(new TargetData(M));
    createStandardFunctionPasses(FPasses, 2);
    FPasses->doInitialization();
    for (Module::iterator I = M->begin(), E = M->end(); I != E; ++I)
      FPasses->run(*I);

    PM.add(new TargetData(M));
    createStandardModulePasses(&PM, 2,
       /* OptimizeSize     */  false,
       /* UnitAtATime      */  true,
       /* UnrollLoops      */  true,
       /* SimplifyLibCalls */  true,
       /* HaveExceptions   */  false,
       /* InliningPass     */  createFunctionInliningPass());

    createStandardLTOPasses(&PM,
      /* Internalize     */ true,
      /* RunInliner      */ true,
      /* VerifyEach      */ false);

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
        InlineFunctionInfo ifi;
        InlineFunction(CI, ifi);
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
            InlineFunctionInfo ifi;
            if ( !InlineFunction(CI, ifi) ){
                std::cerr << __FUNCTION__ << ": Function not inlined\n";
                return false;
            }
            return true;
        }
    }

    std::cout << __FUNCTION__ << ": Function not found\n";
    return false;
}

// CloneFunction wrapper
Function *doCloneFunction(Module *M, Function *Fn, const char *newFnName)
{
    Function *newFn;

    newFn = CloneFunction(Fn);
    newFn->setName(newFnName);
    M->getFunctionList().push_back(newFn);

    return newFn;
}

// CloneFunction wrapper
Function *doCloneFunction(Module *M, const char *FnName, const char *newFnName)
{
    Function *Fn;

    Fn = M->getFunction(FnName);
    assert(Fn);
    return doCloneFunction(M, Fn, newFnName);

}

// Get the first call instruction of Fn in CallerFn
static CallInst *getCallInFn(Function *Fn, Function *CallerFn)
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

// clone and replace the call of a function, inside another function (caller)
Function *CloneAndReplaceFn(Module *M,
                            Function *CallerFn,
                            Function *ReplaceFn,
                            const char *Name,
                            const char *newFnName)
{
    Function *Fn, *newFn;

    // get Function
    Fn = M->getFunction(Name);
    assert(Fn);

    //printf("Uses of %s: %u\n", Fn->getNameStart(), Fn->getNumUses());
    // Clone Caller Function and put it into module
    newFn = doCloneFunction(M, CallerFn, newFnName);

    // get call instruction to replace in new function
    //printf("Uses of %s: %u\n", Fn->getNameStart(), Fn->getNumUses());
    CallInst *CallInstr = getCallInFn(Fn, newFn);
    CallInstr->setOperand(0, ReplaceFn);
    //printf("Uses of %s: %u\n", Fn->getNameStart(), Fn->getNumUses());
    verifyModule(*M, AbortProcessAction, 0);

    return newFn;
}

// wrapper for returning a JIT execution engine
ExecutionEngine *mkJIT(Module *M)
{
    ExecutionEngine *JIT;
    std::string Error;

    JIT = ExecutionEngine::createJIT(M, &Error);
    if (!JIT){
        std::cerr << "ExectionEngine::createJIT:" << Error << "\n";
        exit(1);
    }
    //std::cout << JIT->getTargetData()->getStringRepresentation() << "\n";

    return JIT;
}

/* check that the function is void fn(void) */
static inline int is_hook_type_ok(const Type *Ty)
{
    const FunctionType *FnTy = cast<FunctionType>(Ty);
    return (FnTy->getNumParams() == 0)
           && (FnTy->getReturnType()->isVoidTy());
}

/* find a hook function, and verify its type and that is used once */
static Instruction *get_hook(Module *M, const char *name)
{
    Function *Fn;
    const Type *FnTy;
    Instruction *Hook;

    // get function
    Fn = M->getFunction(name);

    // verify type
    FnTy = Fn->getFunctionType();
    if (!is_hook_type_ok(FnTy)){
        std::cout << "Error in Function Type\n";
        return NULL;
    }

    // get uses, and verify hook is used only once
    Value::use_iterator I = Fn->use_begin(), E = Fn->use_end();
    if (I == E){
        std::cout << "No uses of the hook\n";
        return NULL;
    }
    Hook = cast<Instruction>(*I);
    if (++I != E) {
        std::cout << "More than one uses of the hook\n";
        return NULL;
    }

    return Hook;
}

// similar to above, +verify that hook is used inside Parent
static Instruction *get_hook(Module *M, const char *name, Function *Parent)
{
    Function *Fn;
    Instruction *Hook;
    unsigned int cnt;

    Fn = M->getFunction(name);
    assert(Fn);

    const Type *FnTy = Fn->getFunctionType();
    if (!is_hook_type_ok(FnTy)){
        std::cout << "Error in Function Type\n";
        return NULL;
    }

    Value::use_iterator I = Fn->use_begin();
    Hook = NULL;
    cnt = 0;
    for (; I != Fn->use_end(); ++I){
        Instruction *Ins = cast<Instruction>(*I);
        if (Ins->getParent()->getParent() == Parent){
            cnt++;
            Hook = Ins;
        }
    }

    assert(cnt > 0 && "No uses of Hook found\n");
    assert(cnt < 2 && "More than one uses of Hook found\n");

    return Hook;
}

// core function for replacing a hook with thow basic-blocks
static BasicBlock *__llvm_hook_newbb(Module *M, Instruction *Hook, BasicBlock **BBnext)
{
    BasicBlock *BB = Hook->getParent();
    BasicBlock::iterator BI(Hook);
    BI = BB->getInstList().erase(BI);

    BasicBlock *BB0 = BB->splitBasicBlock(BI, "split_start");
    BasicBlock *BB1 = BB0->splitBasicBlock(BB0->begin(), "split_end");

    // BB0 should have only one (uncond) branch instruction. remove it
    BB0->getInstList().erase(BB0->begin());
    assert(BB0->begin() == BB0->end());

    *BBnext = BB1;
    return BB0;
}


BasicBlock *llvm_hook_newbb(Module *M, const char *name, BasicBlock **BBnext)
{
    Instruction *Hook;

    Hook = get_hook(M,name);
    return __llvm_hook_newbb(M, Hook, BBnext);
}

BasicBlock *llvm_hook_newbb(Module *M, const char *name, Function *Parent,
                            BasicBlock **BBnext)
{
    Instruction *Hook;

    Hook = get_hook(M,name, Parent);
    return __llvm_hook_newbb(M, Hook, BBnext);
}


// update the map with annotations from Module M
// If F is !NULL, then consider only annotations inside F
void Annotations::update(Module *M, Function *Parent)
{
    Function *Fn;
    Value::use_iterator ui;

	// find the annotation function
    Fn = M->getFunction("llvm.var.annotation");
    if (!Fn){
        std::cerr << "no annotations found" << std::endl;
        return;
    }

	// initialize map, if needed
    if (this->map == NULL)
        this->map = new AnnotMap();

	// iterate over uses of the annotation function
    for (ui = Fn->use_begin(); ui != Fn->use_end(); ui++){
        CallInst *CI;
        Value *val;
        GlobalVariable *ptr;
        std::string str;

        // jump trough C++'s hoops
        CI = cast<CallInst>(*ui);
        // if a Parent is included, filter instructions
        if (Parent && CI->getParent()->getParent() != Parent)
            continue;

		// update map
        val = CI->getOperand(1)->stripPointerCasts();
        ptr = cast<GlobalVariable>(CI->getOperand(2)->stripPointerCasts());
        str = cast<ConstantArray>(ptr->getInitializer())->getAsString();
        // insert (the C string)
        (*this->map)[str.c_str()] = val;
    }

	// iterate in all annotation uses, and remove iterations
	// continue, until no annotation is left
    // please, pretty please, fix-me!
    for (;;){
        bool again = false;
        for (ui = Fn->use_begin(); ui != Fn->use_end(); ++ui){
            CallInst *CI;
            CI = cast<CallInst>(*ui);
            if (Parent && CI->getParent()->getParent() != Parent)
                continue;
            CI->eraseFromParent();
            again = true;
            break;
        }
		// no more annotations, nothing else to do
        if (!again)
            break;
    }
}

// get a value from the map
Value *Annotations::getValue(const char *annotation)
{
    AnnotMap::iterator am = this->map->find(annotation);
	if (am == this->map->end())
		return NULL;
    return am->getValue();
}


// dump map
void Annotations::dump()
{
    AnnotMap::iterator am = this->map->begin();

    for (; am != this->map->end(); ++am){
        std::cout << am->getKeyData() << "\n";
    }
}


