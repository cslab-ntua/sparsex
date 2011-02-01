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
#include "llvm/Support/StandardPasses.h"
#include "llvm/LLVMContext.h"

#include "llvm_jit_help.h"

using namespace llvm;

/////////// Utility Functions ///////////

Module *ModuleFromFile(const char *file, LLVMContext Ctx)
{
    std::string Error;
    MemoryBuffer *MB = MemoryBuffer::getFile(file, &Error);
    if (!MB){
        std::cerr << " MemoryBuffer::getFile " << Error << " " << file << "\n";
        exit(1);
    }
    Module *M = ParseBitcodeFile(MB, Ctx, &Error);
    if (!M){
        std::cerr << "ParseBitCodeFile:" << Error << "\n";
        exit(1);
    }

    return M;
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

Function *doCloneFunction(Module *M, Function *Fn, const char *newFnName)
{
    Function *newFn;

    newFn = CloneFunction(Fn);
    newFn->setName(newFnName);
    M->getFunctionList().push_back(newFn);

    return newFn;
}

Function *doCloneFunction(Module *M, const char *FnName, const char *newFnName)
{
    Function *Fn;

    Fn = M->getFunction(FnName);
    assert(Fn);
    return doCloneFunction(M, Fn, newFnName);

}

/*
void PrintAllocas(Value *V)
{
    Value::use_iterator I;
    for (I = V->use_begin(); I != V->use_end(); ++I){
        I->dump();
        AllocaInst *AI = dyn_cast<AllocaInst>(I);
        if (AI){
            std::cout  << "Alloca instruction\n";
        }
    }
}
*/

/* check that the function is void fn(void) */
static inline int is_hook_type_ok(const Type *Ty)
{
    const FunctionType *FnTy = cast<FunctionType>(Ty);
    return (FnTy->getNumParams() == 0)
           && (FnTy->getReturnType()->isVoidTy());
}


static Instruction *get_hook(Module *M, const char *name)
{
    Function *Fn = M->getFunction(name);

    const Type *FnTy = Fn->getFunctionType();
    if (!is_hook_type_ok(FnTy)){
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
    Function *HookFn, *newFn;

    // get Hook Function
    HookFn = M->getFunction(HookName);
    assert(HookFn);

    //printf("Uses of %s: %u\n", HookFn->getNameStart(), HookFn->getNumUses());
    // Clone Caller Function and put it into module

    newFn = doCloneFunction(M, CallerFn, newFnName);

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

static BasicBlock *__llvm_hook_newbb(Module *M, Instruction *Hook, BasicBlock **BBnext)
{
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


// This function replaces a hook function with two basic blocks,
// for the purpose of replacing the hook with other instructions.
// The return basic  block is where the new instructions should be placed,
// while BBnext, is where control flow should continue
BasicBlock *llvm_hook_newbb(Module *M, const char *name, BasicBlock **BBnext)
{
    Instruction *Hook;

    Hook = get_hook(M,name);
    return __llvm_hook_newbb(M, Hook, BBnext);
}

// same with before only hook should be insinde Parent
BasicBlock *llvm_hook_newbb(Module *M, const char *name, Function *Parent,
                            BasicBlock **BBnext)
{
    Instruction *Hook;

    Hook = get_hook(M,name, Parent);
    return __llvm_hook_newbb(M, Hook, BBnext);
}


// update map, with annotations from Module M
// If F is !NULL, then consider only annotations inside F
void Annotations::update(Module *M, Function *Parent)
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
        // if a Parent is included, filter instructions
        if (Parent && CI->getParent()->getParent() != Parent)
            continue;
        val = CI->getOperand(1)->stripPointerCasts();
        ptr = cast<GlobalVariable>(CI->getOperand(2)->stripPointerCasts());
        str = cast<ConstantArray>(ptr->getInitializer())->getAsString();
        // insert (the C string)
        (*this->map)[str.c_str()] = val;
    }

    // please, pretty please fix-me
    for (;;){
        bool erased = false;
        for (ui = Fn->use_begin(); ui != Fn->use_end(); ++ui){
            CallInst *CI;
            CI = cast<CallInst>(*ui);
            if (Parent && CI->getParent()->getParent() != Parent)
                continue;
            CI->eraseFromParent();
            erased = true;
            break;
        }
        if (!erased)
            break;
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

SingleModule::ModMap SingleModule::modules;
SingleModule::RefMap SingleModule::refs;
SingleModule::JitMap SingleModule::jits;

Module *SingleModule::getM(const char *mod_name)
{
    Module *ret;
    ModMap::iterator it;

    it = modules.find(mod_name);
    if (it == modules.end()){
        ret = ModuleFromFile(mod_name);
        modules.insert(std::make_pair(mod_name, ret));
        refs.insert(std::make_pair(ret, 1));
    } else {
        ret = it->second;
        refs[ret] += 1;
    }

    return ret;
}

ExecutionEngine *mkJIT(Module *M)
{
    ExecutionEngine *JIT;
    std::string Error;

    JIT = ExecutionEngine::createJIT(M, &Error);
    if (!JIT){
        std::cerr << "ExectionEngine::createJIT:" << Error << "\n";
        exit(1);
    }

    return JIT;
}


ExecutionEngine *SingleModule::getJIT(Module *M)
{
    ExecutionEngine *ret;
    JitMap::iterator it;

    it = jits.find(M);
    if (it == jits.end()){
        ret = mkJIT(M);
        jits.insert(std::make_pair(M, ret));
    } else {
        ret = it->second;
    }

    return ret;
}


#if 0
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

extern JitSymbolTable *__jitSymbolTable;

unsigned long getFnSize(const char *FnName)
{
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

    return 0;
}

void printSymTable()
{
    unsigned int i;
    for (i=0; i<__jitSymbolTable->NumSymbols; i++){
        JitSymbolEntry *entry = &__jitSymbolTable->Symbols[i];
        std::cout << entry->FnName << " : " << entry->FnStart << " " << entry->FnSize << "\n";
    }
}

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
