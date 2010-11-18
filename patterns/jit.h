#ifndef CSX_JIT_H__
#define CSX_JIT_H__

#include "csx.h"

#include "llvm/Type.h"
#include "llvm/Module.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Function.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm_jit_help.h"

using namespace llvm;

namespace csx {

class CsxJit {
public:
    CsxManager *CsxMg;
    Module *M;
    IRBuilder<> *Bld;
    //ModuleProvider *MP;
    //ExecutionEngine *JIT;

    Function *UlGet;
    Function *DecodeF;
    Function *SpmvF;
    Function *FailF;
    Function *PrintYXV;
    Function *AlignF;
    Function *TestBitF;

    Value *Zero8;
    Value *Zero32;
    Value *Zero64;
    Value *One8;
    Value *One32;
    Value *One64;
    Value *Three64;

    Value *YrPtr;
    Value *MyxPtr;
    Value *Xptr;
    Value *Yptr;
    Value *YindxPtr;
    Value *Vptr;
    Value *CtlPtr;
    Value *SizePtr;
    Value *FlagsPtr;

    Annotations annotations;

    CsxJit(CsxManager *CsxMg, unsigned int thread_id=0);

    void doPrint(Value *Myx=NULL, Value *Yindx=NULL);
    void doIncV();
    void doMul(Value *Myx=NULL, Value *Yindx=NULL);
    void doStoreYr();
    void doOp(Value *MyX=NULL, Value *Yindx=NULL);

    void doNewRowHook();
    void doBodyHook(std::ostream &os);
    void doHooks(std::ostream &os);
    void *doJit();

    void doDeltaAddMyx(int delta_bytes);

    void DeltaCase(BasicBlock *BB,        // case
	           BasicBlock *BB_lentry, // loop entry
	           BasicBlock *BB_lbody,  // loop body
	           BasicBlock *BB_exit,   // final exit
	           int delta_bytes);

    void HorizCase(BasicBlock *BB,
	           BasicBlock *BB_lbody,
	           BasicBlock *BB_lexit,
                   BasicBlock *BB_exit,
                   int delta_size);

    void VertCase(BasicBlock *BB,
                  BasicBlock *BB_lbody,
                  BasicBlock *BB_exit,
                  int delta_size);

    void DiagCase(BasicBlock *BB,
                  BasicBlock *BB_lbody,
                  BasicBlock *BB_exit,
                  int delta_size,
                  bool reversed);

    void BlockRowCaseRolled(BasicBlock *BB,
	                    BasicBlock *BB_lbody,
                            BasicBlock *BB_exit,
                            int r, int c);

    void BlockColCaseRolled(BasicBlock *BB,
	                    BasicBlock *BB_lbody,
                            BasicBlock *BB_exit,
                            int r, int c);

    void BlockRowCaseUnrolled(BasicBlock *BB,
                              BasicBlock *BB_exit,
                              int r, int c);

    void BlockColCaseUnrolled(BasicBlock *BB,
                              BasicBlock *BB_exit,
                              int r, int c);
};



} // end csx namespace

#endif // CSX_JIT_H__
