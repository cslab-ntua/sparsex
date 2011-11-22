/*
 * jit.cc -- Just In Time compilation routines.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2011, Vasileios Karakasis
 * Copyright (C) 2010-2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <iostream>
#include <sstream>
#include <cassert>

#include "llvm/Analysis/Verifier.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/Target/TargetSelect.h"

#include "spm.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"

#include "llvm_jit_help.h"

using namespace llvm;
using namespace csx;

#ifndef CSX_TEMPLATE
#   define CSX_TEMPLATE "csx_llvm_tmpl.llvm.bc"
#endif

///< State shared between all CsxJits.
///< We maintain a single Module and Context for LLVM.
static LLVMContext Ctx_;
static Module *Mod_ = NULL;

void CsxJitInitGlobal(void)
{
    //std::cout << __FUNCTION__ << ": One-time initialization" << "\n";
    InitializeNativeTarget();
    Mod_ = ModuleFromFile(CSX_TEMPLATE, Ctx_);
}

void CsxJitOptmize(void)
{
    assert(Mod_ && "module not initiated");
    doOptimize(Mod_);
}

///< Since LLVM Module is shared, Use a common JIT.
static ExecutionEngine *DoJIT(void)
{
    static ExecutionEngine *ee = NULL;

    if (ee == NULL)
        ee = mkJIT(Mod_);

    return ee;
}

LLVMContext &CsxJit::GetLLVMCtx()
{
    assert(Mod_ && "module not initiated");
    return Ctx_;
}

CsxJit::CsxJit(CsxManager *csx_mg, unsigned int tid) : CsxMg(csx_mg)
{
    assert(Mod_ && "module not initiated");
    this->M = Mod_;
    this->Bld = new IRBuilder<>(GetLLVMCtx());

    ///< Create a new spmv function for the thread.
    std::ostringstream str_stream;
    str_stream << "csx_spmv_" << tid;
    this->SpmvF = doCloneFunction(M, "csx_spmv_template",
                                  str_stream.str().c_str());
    //str_stream << ".llvm.bc";
    //ModuleToFile(this->M, str_stream.str().c_str());

    ///< Load helper functions from the compiled template.
    this->PrintYXV = M->getFunction("print_yxv");
    this->FailF    = M->getFunction("fail");
    this->AlignF   = M->getFunction("align_ptr");
    this->TestBitF = M->getFunction("test_bit");
    this->UlGet    = M->getFunction("ul_get");

    ///< Find annotated values from the template.
    Annotations annotations;
    annotations.update(M, this->SpmvF);
    //this->annotations.dump();
    this->YrPtr    = annotations.getValue("spmv::yr");
    this->MyxPtr   = annotations.getValue("spmv::myx");
    this->Xptr     = annotations.getValue("spmv::x");
    this->Yptr     = annotations.getValue("spmv::y");
    this->YindxPtr = annotations.getValue("spmv::y_indx");
    this->Vptr     = annotations.getValue("spmv::v");
    this->CtlPtr   = annotations.getValue("spmv::ctl");
    this->SizePtr  = annotations.getValue("spmv::size");
    this->FlagsPtr = annotations.getValue("spmv::flags");

    ///< Initialize needed constants.
    this->Zero8   = ConstantInt::get(Type::getInt8Ty (GetLLVMCtx()), 0);
    this->Zero32  = ConstantInt::get(Type::getInt32Ty(GetLLVMCtx()), 0);
    this->Zero64  = ConstantInt::get(Type::getInt64Ty(GetLLVMCtx()), 0);
    this->One8    = ConstantInt::get(Type::getInt8Ty (GetLLVMCtx()), 1);
    this->One32   = ConstantInt::get(Type::getInt32Ty(GetLLVMCtx()), 1);
    this->One64   = ConstantInt::get(Type::getInt64Ty(GetLLVMCtx()), 1);
    this->Three64 = ConstantInt::get(Type::getInt64Ty(GetLLVMCtx()), 3);
}

void CsxJit::DoStoreYr()
{
    Value *Yr, *Yindx, *YiPtr, *Yi;

    Yr = Bld->CreateLoad(YrPtr, "yr");

    YiPtr = Bld->CreateLoad(Yptr, "y");
    Yindx = Bld->CreateLoad(YindxPtr, "y_indx");
    YiPtr = Bld->CreateGEP(YiPtr, Yindx, "yi");

    Yi = Bld->CreateLoad(YiPtr, "yi");
    Yi = Bld->CreateAdd(Yi, Yr);

    Bld->CreateStore(Yi, YiPtr);
}

void CsxJit::DoNewRowHook()
{
    BasicBlock *BB, *BB_next;
    Value *v;

    ///< New row.
    BB = llvm_hook_newbb(M, "__new_row_hook", SpmvF, &BB_next);
    Bld->SetInsertPoint(BB);
    DoStoreYr();

    if (!CsxMg->HasRowJmps()) {
        v = Bld->CreateLoad(YindxPtr, "y_indx");
        v = Bld->CreateAdd(v, One64, "y_indx_inc");
        Bld->CreateStore(v, YindxPtr);
        Bld->CreateBr(BB_next);
    } else {
        BasicBlock *BB_rjmp, *BB_rend;
        Value *RJmpBit, *Yindx, *Test;
        Value *Ul;
        PHINode *YindxAdd;

        BB_rjmp = BasicBlock::Create(GetLLVMCtx(), "rjmp",
                                     BB->getParent(), BB_next);
        BB_rend = BasicBlock::Create(GetLLVMCtx(), "rend",
                                     BB->getParent(), BB_next);
        RJmpBit = ConstantInt::get(Type::getInt32Ty(GetLLVMCtx()),
                                   CTL_RJMP_BIT);
        Yindx = Bld->CreateLoad(YindxPtr, "y_indx");
        Test = Bld->CreateCall2(TestBitF, FlagsPtr, RJmpBit);
        Test = Bld->CreateICmpEQ(Test, Zero32, "bit_test");
        Bld->CreateCondBr(Test, BB_rend, BB_rjmp);
        Bld->SetInsertPoint(BB_rjmp);
        Ul = Bld->CreateCall(UlGet, CtlPtr);
        Bld->CreateBr(BB_rend);

        ///< Common end.
        Bld->SetInsertPoint(BB_rend);
        YindxAdd = Bld->CreatePHI(Type::getInt64Ty(GetLLVMCtx()), "yindx_add");
        YindxAdd->addIncoming(One64, BB);
        YindxAdd->addIncoming(Ul, BB_rjmp);

        v = Bld->CreateAdd(YindxAdd, Yindx);
        Bld->CreateStore(v, YindxPtr);
        Bld->CreateBr(BB_next);
    }
}

void CsxJit::HorizCase(BasicBlock *BB,
                       BasicBlock *BB_lbody, BasicBlock *BB_lexit,
                       BasicBlock *BB_exit,
                       int delta_size)
{
    Value *Size, *Delta, *Myx0, *newMyx, *NextCnt, *Test;
    PHINode *Myx, *Cnt;

    Delta = ConstantInt::get(Type::getInt64Ty(GetLLVMCtx()), delta_size);

    Bld->SetInsertPoint(BB);
    Size = Bld->CreateLoad(SizePtr, "size");
    Myx0 = Bld->CreateLoad(MyxPtr, "myx0");
    Bld->CreateBr(BB_lbody);

    ///< Body
    Bld->SetInsertPoint(BB_lbody);
    Cnt = Bld->CreatePHI(Type::getInt8Ty(GetLLVMCtx()), "cnt");
    Myx = Bld->CreatePHI(Myx0->getType(), "myx");
    DoOp(Myx);

    newMyx = Bld->CreateGEP(Myx, Delta, "new_myx");

    NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");
    Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
    Bld->CreateCondBr(Test, BB_lexit, BB_lbody);

    Cnt->addIncoming(Zero8, BB);
    Cnt->addIncoming(NextCnt, BB_lbody);

    Myx->addIncoming(Myx0, BB);
    Myx->addIncoming(newMyx, BB_lbody);

    ///< Exit
    Bld->SetInsertPoint(BB_lexit);
    Bld->CreateStore(Myx, MyxPtr);
    Bld->CreateBr(BB_exit);
}

void CsxJit::VertCase(BasicBlock *BB,
                      BasicBlock *BB_lbody,
                      BasicBlock *BB_exit,
                      int delta_size)
{
    Value *Size, *Delta, *Yindx0, *YindxAdd, *NextCnt, *Test;
    PHINode *Yindx, *Cnt;

    Delta = ConstantInt::get(Type::getInt64Ty(GetLLVMCtx()), delta_size);

    Bld->SetInsertPoint(BB);
    Size = Bld->CreateLoad(SizePtr, "size");
    Yindx0 = Bld->CreateLoad(YindxPtr);
    Bld->CreateBr(BB_lbody);

    ///< Body
    Bld->SetInsertPoint(BB_lbody);
    Cnt =   Bld->CreatePHI(Type::getInt8Ty (GetLLVMCtx()), "cnt");
    Yindx = Bld->CreatePHI(Type::getInt64Ty(GetLLVMCtx()), "yindx");

    DoOp(NULL, Yindx);

    YindxAdd = Bld->CreateAdd(Yindx, Delta);
    NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");
    Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
    Bld->CreateCondBr(Test, BB_exit, BB_lbody);

    Cnt->addIncoming(Zero8, BB);
    Cnt->addIncoming(NextCnt, BB_lbody);

    Yindx->addIncoming(Yindx0, BB);
    Yindx->addIncoming(YindxAdd, BB_lbody);
}

void CsxJit::DiagCase(BasicBlock *BB,
                      BasicBlock *BB_lbody,
                      BasicBlock *BB_exit,
                      int delta_size,
                      bool reversed)
{
    Value *Size, *D, *minusD, *Test;
    PHINode *Myx, *Yindx, *Cnt;
    Value *Myx0, *Yindx0;
    Value *newMyx, *YindxAdd, *NextCnt;

    D = ConstantInt::get(Type::getInt64Ty(GetLLVMCtx()), delta_size);
    minusD = ConstantInt::get(Type::getInt64Ty(GetLLVMCtx()), -delta_size);

    Bld->SetInsertPoint(BB);
    Size = Bld->CreateLoad(SizePtr, "size");
    Myx0 = Bld->CreateLoad(MyxPtr, "myx0");
    Yindx0 = Bld->CreateLoad(YindxPtr);
    Bld->CreateBr(BB_lbody);

    ///< Body
    Bld->SetInsertPoint(BB_lbody);
    Cnt = Bld->CreatePHI(Type::getInt8Ty(GetLLVMCtx()), "cnt");
    Yindx = Bld->CreatePHI(Type::getInt64Ty(GetLLVMCtx()), "yindx");
    Myx = Bld->CreatePHI(Myx0->getType(), "myx");

    DoOp(Myx, Yindx);

    YindxAdd = Bld->CreateAdd(Yindx, D);
    newMyx = reversed ? Bld->CreateGEP(Myx, minusD) : Bld->CreateGEP(Myx, D);
    NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");
    Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
    Bld->CreateCondBr(Test, BB_exit, BB_lbody);

    Cnt->addIncoming(Zero8, BB);
    Cnt->addIncoming(NextCnt, BB_lbody);

    Myx->addIncoming(Myx0, BB);
    Myx->addIncoming(newMyx, BB_lbody);

    Yindx->addIncoming(Yindx0, BB);
    Yindx->addIncoming(YindxAdd, BB_lbody);
}

void CsxJit::BlockRowCaseRolled(BasicBlock *BB,
                                BasicBlock *BB_lbody,
                                BasicBlock *BB_exit,
                                int r, int c)
{
    Value *Myx0, *NewMyx, *NextCnt, *Size, *Test, *Yindx;
    PHINode *Myx, *Cnt;

    ///< Initialization
    Bld->SetInsertPoint(BB);
    Size = ConstantInt::get(Type::getInt8Ty(GetLLVMCtx()), c);
    Myx0 = Bld->CreateLoad(MyxPtr);
    Bld->CreateBr(BB_lbody);

    ///< Body
    Bld->SetInsertPoint(BB_lbody);
    Cnt = Bld->CreatePHI(Type::getInt8Ty(GetLLVMCtx()), "cnt");
    Myx = Bld->CreatePHI(Myx0->getType(), "myx");
    Yindx = Bld->CreateLoad(YindxPtr, "yindx");
    for (int i = 0; i < r; i++) {
        DoOp(Myx, Yindx);
        Yindx = Bld->CreateAdd(Yindx, One64);
    }

    NewMyx = Bld->CreateGEP(Myx, One64);
    NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");
    Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
    Bld->CreateCondBr(Test, BB_exit, BB_lbody);

    Cnt->addIncoming(Zero8, BB);
    Cnt->addIncoming(NextCnt, BB_lbody);

    Myx->addIncoming(Myx0, BB);
    Myx->addIncoming(NewMyx, BB_lbody);
}

void CsxJit::BlockColCaseRolled(BasicBlock *BB,
                                BasicBlock *BB_lbody,
                                BasicBlock *BB_exit,
                                int r, int c)
{
    Value *Yindx0, *NewYindx, *NextCnt, *Size, *Test, *Myx;
    PHINode *Yindx, *Cnt;

    ///< Initialization
    Bld->SetInsertPoint(BB);
    Size = ConstantInt::get(Type::getInt8Ty(GetLLVMCtx()), r);
    Yindx0 = Bld->CreateLoad(YindxPtr);
    Bld->CreateBr(BB_lbody);

    ///< Body
    Bld->SetInsertPoint(BB_lbody);
    Cnt = Bld->CreatePHI(Type::getInt8Ty(GetLLVMCtx()), "cnt");
    Yindx = Bld->CreatePHI(Yindx0->getType(), "yindx");
    Myx = Bld->CreateLoad(MyxPtr,"myx");
    for (int i = 0; i < c; i++) {
        DoOp(Myx, Yindx);
        Myx = Bld->CreateGEP(Myx, One64);
    }

    NewYindx = Bld->CreateAdd(Yindx, One64);
    NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");
    Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
    Bld->CreateCondBr(Test, BB_exit, BB_lbody);

    Cnt->addIncoming(Zero8, BB);
    Cnt->addIncoming(NextCnt, BB_lbody);

    Yindx->addIncoming(Yindx0, BB);
    Yindx->addIncoming(NewYindx, BB_lbody);
}

void CsxJit::BlockRowCaseUnrolled(BasicBlock *BB,
                                  BasicBlock *BB_exit,
                                  int r, int c)
{
    Value **Myx = new Value*[c+1];
    Value **Yindx = new Value*[r+1];

    Bld->SetInsertPoint(BB);
    Myx[0] = Bld->CreateLoad(MyxPtr);

#ifdef NO_EXTRA_LOADS
    Yindx[0] = Bld->CreateLoad(YindxPtr);
    for (int i = 0; i < c; i++)
        Myx[i+1] = Bld->CreateGEP(Myx[i], One64);

    for (int j = 0; j < r; j++)
        Yindx[j+1] = Bld->CreateAdd(Yindx[j], One64);

    for (int i = 0; i < c; i++)
        for (int j = 0; j < r; j++)
            DoOp(Myx[i], Yindx[j]);
#else
    ///< Elements in block-row types are stored column-wise
    for (int i = 0; i < c; i++) {
        Yindx[0] = Bld->CreateLoad(YindxPtr);
        for (int j = 0; j < r; j++) {
            DoOp(Myx[i], Yindx[j]);
            Yindx[j+1] = Bld->CreateAdd(Yindx[j], One64);
        }
        Myx[i+1] = Bld->CreateGEP(Myx[i], One64);
    }
#endif
    Bld->CreateBr(BB_exit);
}

void CsxJit::BlockColCaseUnrolled(BasicBlock *BB,
                                  BasicBlock *BB_exit,
                                  int r, int c)
{
    Value **Myx = new Value*[c+1];
    Value **Yindx = new Value*[r+1];

    Bld->SetInsertPoint(BB);
    ///< Elements in block-col types are stored row-wise
    Yindx[0] = Bld->CreateLoad(YindxPtr);
#ifdef NO_EXTRA_LOADS
    Myx[0] = Bld->CreateLoad(MyxPtr);
    for (int j = 0; j < c; j++)
        Myx[j+1] = Bld->CreateGEP(Myx[j], One64);

    for (int i = 0; i < r; i++)
        Yindx[i+1] = Bld->CreateAdd(Yindx[i], One64);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            DoOp(Myx[j], Yindx[i]);
#else
    for (int i = 0; i < r; i++) {
        Myx[0] = Bld->CreateLoad(MyxPtr);
        for (int j = 0; j < c; j++) {
            DoOp(Myx[j], Yindx[i]);
            Myx[j+1] = Bld->CreateGEP(Myx[j], One64);
        }
        Yindx[i+1] = Bld->CreateAdd(Yindx[i], One64);
    }
#endif

    Bld->CreateBr(BB_exit);
}

void CsxJit::DoPrint(Value *Myx, Value *Yindx)
{
    Value *X, *Xindx;
    if (Myx == NULL)
        Myx = Bld->CreateLoad(MyxPtr);

    Myx = Bld->CreatePtrToInt(Myx, Type::getInt64Ty(GetLLVMCtx()), "myx_int");
    X = Bld->CreateLoad(Xptr);
    X = Bld->CreatePtrToInt(X, Type::getInt64Ty(GetLLVMCtx()), "x_int");
    Xindx = Bld->CreateSub(Myx, X);
    Xindx = Bld->CreateAShr(Xindx, Three64);

    if (Yindx == NULL)
        Yindx = Bld->CreateLoad(YindxPtr);

    Value *V;
    V = Bld->CreateLoad(Bld->CreateLoad(Vptr));

    Bld->CreateCall3(PrintYXV, Yindx, Xindx, V);
}

void CsxJit::DoMul(Value *Myx, Value *Yindx)
{
    Value *V, *X;

    if (Myx == NULL)
        Myx = Bld->CreateLoad(MyxPtr);

    X = Bld->CreateLoad(Myx, "x");
    V = Bld->CreateLoad(Bld->CreateLoad(Vptr, "v_ptr"), "val");
    V = Bld->CreateMul(V, X, "mul");

    if (Yindx == NULL) {
        Value *Yr;
        ///< Use Yr to store the result.
        Yr = Bld->CreateLoad(YrPtr, "yr");
        Yr = Bld->CreateAdd(Yr, V, "yr_add");
        Bld->CreateStore(Yr, YrPtr);
    } else {
        Value *Yi, *YiPtr;
        YiPtr = Bld->CreateLoad(Yptr, "y");
        YiPtr = Bld->CreateGEP(YiPtr, Yindx);
        Yi = Bld->CreateLoad(YiPtr, "yi");
        Yi = Bld->CreateAdd(Yi, V, "new_yi");
        Bld->CreateStore(Yi, YiPtr);
    }
}

void CsxJit::DoIncV()
{
    Value *V = Bld->CreateLoad(Vptr);
    Value *newV = Bld->CreateGEP(V, One64);
    Bld->CreateStore(newV, Vptr);
}

void CsxJit::DoOp(Value *Myx, Value *Yindx)
{
    DoMul(Myx, Yindx);
    DoIncV();
}

void CsxJit::DoDeltaAddMyx(int delta_bytes)
{
    Function *F = NULL;

    switch (delta_bytes) {
        case 1:
            F = M->getFunction("u8_get");
            break;

        case 2:
            F = M->getFunction("u16_get");
            break;

        case 4:
            F = M->getFunction("u32_get");
            break;

        case 8:
            F = M->getFunction("u64_get");
            break;

        default:
            assert(false);
    }
    assert(F && "Failed to find uxx_get() function");

    Value *Myx = Bld->CreateLoad(MyxPtr, "myx");
    Value *MyxAdd = Bld->CreateCall(F, CtlPtr, "myx_add");
    Value *newMyx = Bld->CreateGEP(Myx, MyxAdd, "newmyx");
    Bld->CreateStore(newMyx, MyxPtr);
}

void CsxJit::DeltaCase(BasicBlock *BB,
                       BasicBlock *BB_entry, BasicBlock *BB_body,
                       BasicBlock *BB_exit,
                       int delta_bytes)
{
    Value *Align, *Size, *Test, *NextCnt;
    PHINode *Cnt;

    Bld->SetInsertPoint(BB);
    ///< Align ctl.
    if (delta_bytes > 1) {
        Align = ConstantInt::get(Type::getInt32Ty(GetLLVMCtx()),delta_bytes);
        Bld->CreateCall2(AlignF, CtlPtr, Align);
    }
    Size = Bld->CreateLoad(SizePtr, "size");
    Bld->CreateBr(BB_entry);

    ///< Entry
    Bld->SetInsertPoint(BB_entry);

    DoOp();

    Test = Bld->CreateICmpUGT(Size, One8);
    Bld->CreateCondBr(Test, BB_body, BB_exit);

    ///< Body
    Bld->SetInsertPoint(BB_body);
    Cnt = Bld->CreatePHI(Type::getInt8Ty(GetLLVMCtx()), "cnt");
    DoDeltaAddMyx(delta_bytes);
    NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");

    DoOp();

    Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
    Bld->CreateCondBr(Test, BB_exit, BB_body);

    Cnt->addIncoming(One8, BB_entry);
    Cnt->addIncoming(NextCnt, BB_body);
}

void CsxJit::DoBodyHook(std::ostream &os)
{
    BasicBlock *BB, *BB_next, *BB_default, *BB_case;
    Value *PatternMask;
    Value *v;
    uint64_t delta;
    SpmIterOrder type;

    BB = llvm_hook_newbb(M, "__body_hook", SpmvF, &BB_next);

    ///< Get pattern for switch instruction.
    Bld->SetInsertPoint(BB);
    PatternMask = ConstantInt::get(Type::getInt8Ty(GetLLVMCtx()),
                                   CTL_PATTERN_MASK);
    v = Bld->CreateLoad(FlagsPtr, "flags");
    v = Bld->CreateAnd(PatternMask, v, "pattern");

    ///< Switch default block (call the fail function).
    BB_default = BasicBlock::Create(GetLLVMCtx(), "default", BB->getParent(),
                                    BB_next);
    Bld->SetInsertPoint(BB_default);
    Bld->CreateCall(FailF);
    Bld->CreateBr(BB_next);

    ///< Switch instruction.
    SwitchInst *Switch;
    Bld->SetInsertPoint(BB);
    Switch = Bld->CreateSwitch(v, BB_default, CsxMg->patterns.size());

    ///< Fill up switch, by iterating given patterns.
    CsxManager::PatMap::iterator pat_i = CsxMg->patterns.begin();
    BasicBlock *BB_lentry, *BB_lbody, *BB_lexit;
    for ( ; pat_i !=  CsxMg->patterns.end(); ++pat_i) {
        ///< Alocate case + loop BBs.
        BB_case = BasicBlock::Create(GetLLVMCtx(), "case", BB->getParent(),
                                     BB_default);
        type  = static_cast<SpmIterOrder>(pat_i->first / CSX_PID_OFFSET);
        delta = pat_i->first % CSX_PID_OFFSET;
        switch (type) {
        ///< Deltas
        case 0:
            assert(delta ==  8 ||
                   delta == 16 ||
                   delta == 32 ||
                   delta == 64);
            os << "type:DELTA size:" << delta << " nnz:"
               << pat_i->second.nr << std::endl;
            BB_lentry = BasicBlock::Create(GetLLVMCtx(), "lentry",
                                           BB->getParent(), BB_default);
            BB_lbody  = BasicBlock::Create(GetLLVMCtx(), "lbody",
                                           BB->getParent(), BB_default);
            DeltaCase(BB_case,
                      BB_lentry, BB_lbody,
                      BB_next,
                      delta / 8);
            break;

        ///< Horizontal
        case HORIZONTAL:
            os << "type:HORIZONTAL delta:" << delta
               << " nnz:" << pat_i->second.nr << std::endl;
            BB_lbody = BasicBlock::Create(GetLLVMCtx(), "lbody",
                                          BB->getParent(), BB_default);
            BB_lexit = BasicBlock::Create(GetLLVMCtx(), "lexit",
                                          BB->getParent(), BB_default);
            HorizCase(BB_case,
                      BB_lbody, BB_lexit,
                      BB_next,
                      delta);
            break;

        ///< Vertical
        case VERTICAL:
            os << "type:VERTICAL delta:" << delta
               << " nnz:" << pat_i->second.nr << std::endl;
            BB_lbody = BasicBlock::Create(GetLLVMCtx(), "lbody",
                                          BB->getParent(), BB_default);
            VertCase(BB_case,
                     BB_lbody,
                     BB_next,
                     delta);
            break;

        ///< Diagonal
        case DIAGONAL:
            os << "type:DIAGONAL delta:" << delta
               << " nnz:" << pat_i->second.nr << std::endl;
            BB_lbody = BasicBlock::Create(GetLLVMCtx(), "lbody",
                                          BB->getParent(), BB_default);
            DiagCase(BB_case,
                     BB_lbody,
                     BB_next,
                     delta,
                     false);
            break;

        ///< Reverse Diagonal
        case REV_DIAGONAL:
            os << "type:REV_DIAGONAL delta:" << delta
               << " nnz:" << pat_i->second.nr << std::endl;
            BB_lbody = BasicBlock::Create(GetLLVMCtx(), "lbody",
                                          BB->getParent(), BB_default);
            DiagCase(BB_case,
                     BB_lbody,
                     BB_next,
                     delta,
                     true);
            break;

        ///< Row blocks
        case BLOCK_TYPE_START ... BLOCK_COL_START - 1:
            os << "type:" << SpmTypesNames[type]
               << " dim:" << (type - BLOCK_TYPE_START) << "x" << delta 
               << " nnz:" << pat_i->second.nr << std::endl;
            BB_lbody = BasicBlock::Create(GetLLVMCtx(), "lbody",
                                          BB->getParent(), BB_default);
            BlockRowCaseRolled(BB_case, BB_lbody, BB_next,
                               type - BLOCK_TYPE_START, delta);
            break;
            
        ///< Column Blocks
        case BLOCK_COL_START ... BLOCK_TYPE_END:
            os << "type:" << SpmTypesNames[type] 
               << " dim:" << delta << "x" << (type - BLOCK_COL_START)
               << " nnz:" <<  pat_i->second.nr << std::endl;
            BB_lbody = BasicBlock::Create(GetLLVMCtx(), "lbody",
                                          BB->getParent(), BB_default);
            BlockColCaseRolled(BB_case, BB_lbody, BB_next, delta,
                               type - BLOCK_COL_START);
            /*BlockColCaseUnrolled(BB_case, BB_next, delta,
                                   type - BLOCK_COL_START);*/
            break;
            
        default:
            assert(false && "there is no such type");
        }

        Switch->addCase(
            ConstantInt::get(Type::getInt8Ty(GetLLVMCtx()), pat_i->second.flag),
            BB_case
        );
    }
}

void CsxJit::GenCode(std::ostream &os)
{
    DoNewRowHook();
    DoBodyHook(os);
}

void *CsxJit::GetSpmvFn()
{
    ExecutionEngine *JIT;

    verifyModule(*M, AbortProcessAction, 0);
    //ModuleToFile(M, "M.llvm.bc");
    //DoOptimize(M);
    //M->dump();
    //std::cerr << "Generating Function\n";
    JIT = DoJIT(); //SingleModule::getJIT(M);
    return JIT->getPointerToFunction(SpmvF);
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
