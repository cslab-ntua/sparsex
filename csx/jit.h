#ifndef CSX_JIT_H__
#define CSX_JIT_H__

#include "csx.h"

#include "llvm/Type.h"
#include "llvm/Module.h"
#include "llvm/LLVMContext.h"
#include "llvm/Support/IRBuilder.h"
#include "llvm/Function.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm_jit_help.h"

using namespace llvm;

namespace csx {

// CsxJit is responsible for generating a spmv version for the current CSX
// matrix. Because we don't want to generate all the code from scratch using
// LLVM calls, we use a template written in C, and compile it into LLVM
// bytecode.
//
// The template is stored in CSX_TEMPLATE and contains:
//  - various helper functions
//    (e.g., print_yxv(), fail(), align_ptr(), ...)
//  - a template function called csx_spmv_template.
//
// We use the template function to generate our code:
//  - we clone it, and create a separate function for each thread
//  - we use annotations to mark specific values that we need to access
//  inside the function.
//  - The function contains hooks that are overwritten, with matrix-specific
//  generated code
//
// see also: csx_llvm_tmpl.c

class CsxJit {
public:

    // main state
    LLVMContext Ctx;
    CsxManager  *CsxMg;
    Module      *M;
    IRBuilder<> *Bld;
    //ModuleProvider *MP;
    //ExecutionEngine *JIT;

    // Helper functions that are loaded from the template
    Function
      *UlGet,    // extract a variable-length unsigned long from encoded data
      *FailF,    // a function to indicate that something went wrong
      *PrintYXV, // print a triplet of Y, X, Value
      *AlignF,   // align a pointer to a specific boundary
      *TestBitF, // test if a bit is set
      *SpmvF;    // a copied instance of the template's spmv function

    // Commonly used Constants
    Value *Zero8, *Zero32, *Zero64, *One8, *One32, *One64, *Three64;

    // Values for handling the template function's local variables
    // Note: We load and store into the variables stack addresses,
    //       We leave it to optimization to promote mem to regsiters.
    Value
       *Xptr,      // X:     X array
       *Yptr,      // Y:     Y array
       *Vptr,      // V:     V array
       *YrPtr,     // Yr:    current partial value of y
       *CtlPtr,    // Ctl:   Ctl array (encoded index information)
       *MyxPtr,    // MyX:   current position in X (pointer)
       *YindxPtr,  // Yindx: Y index
       *SizePtr,   // Size:  size of current unit (taken from Ctl)
       *FlagsPtr;  // Flags: flags of current unit (taken from Ctl)


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

    void DeltaCase(BasicBlock *BB,    // case
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

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
