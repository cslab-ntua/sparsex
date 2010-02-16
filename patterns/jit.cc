#include <iostream>
#include <sstream>
#include <cassert>

#include "llvm/Analysis/Verifier.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ModuleProvider.h"

#include "spm.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"

#include "llvm_jit_help.h"

using namespace llvm;
using namespace csx;

CsxJit::CsxJit(CsxManager *_ctl_mg, unsigned int tid) : CsxMg(_ctl_mg)
{
	std::ostringstream str_stream;

	this->M = SingleModule::getM("csx_llvm_tmpl.llvm.bc");
	this->Bld = new IRBuilder<>();

	str_stream << "csx_spmv_" << tid;
	this->SpmvF = doCloneFunction(M, "csx_spmv_template", str_stream.str().c_str());

	this->annotations.update(M, this->SpmvF);
	//this->annotations.dump();

	//str_stream << ".llvm.bc";
	//ModuleToFile(this->M, str_stream.str().c_str());

	this->PrintYXV = M->getFunction("print_yxv");
	this->FailF = M->getFunction("fail");
	this->AlignF = M->getFunction("align_ptr");
	this->TestBitF = M->getFunction("test_bit");
	this->UlGet = M->getFunction("ul_get");

	this->YrPtr = annotations.getValue("spmv::yr");
	this->MyxPtr = annotations.getValue("spmv::myx");
	this->Xptr = annotations.getValue("spmv::x");
	this->Yptr = annotations.getValue("spmv::y");
	this->YindxPtr = annotations.getValue("spmv::y_indx");
	this->Vptr = annotations.getValue("spmv::v");
	this->CtlPtr = annotations.getValue("spmv::ctl");
	this->SizePtr = annotations.getValue("spmv::size");
	this->FlagsPtr = annotations.getValue("spmv::flags");

	this->Zero8 = ConstantInt::get(Type::Int8Ty, 0);
	this->Zero32 = ConstantInt::get(Type::Int32Ty, 0);
	this->Zero64 = ConstantInt::get(Type::Int64Ty, 0);
	this->One8 = ConstantInt::get(Type::Int8Ty, 1);
	this->One32 = ConstantInt::get(Type::Int32Ty, 1);
	this->One64 = ConstantInt::get(Type::Int64Ty, 1);
	this->Three64 = ConstantInt::get(Type::Int64Ty, 3);
}

void CsxJit::doStoreYr()
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

void CsxJit::doNewRowHook()
{
	BasicBlock *BB, *BB_next;
	Value *v;

	// new row
	BB = llvm_hook_newbb(M, "__new_row_hook", SpmvF, &BB_next);

	Bld->SetInsertPoint(BB);
	doStoreYr();

	if (!CsxMg->row_jmps){
		v = Bld->CreateLoad(YindxPtr, "y_indx");
		v = Bld->CreateAdd(v, One64, "y_indx_inc");
		Bld->CreateStore(v, YindxPtr);
		Bld->CreateBr(BB_next);
	} else {
		BasicBlock *BB_rjmp, *BB_rend;
		Value *RJmpBit, *Yindx, *Test;
		Value *Ul;
		PHINode *YindxAdd;

		BB_rjmp = BasicBlock::Create("rjmp", BB->getParent(), BB_next);
		BB_rend = BasicBlock::Create("rend", BB->getParent(), BB_next);
		RJmpBit = ConstantInt::get(Type::Int32Ty, CTL_RJMP_BIT);

		Yindx = Bld->CreateLoad(YindxPtr, "y_indx");
		Test = Bld->CreateCall2(TestBitF, FlagsPtr, RJmpBit);
		Test = Bld->CreateICmpEQ(Test, Zero32, "bit_test");
		Bld->CreateCondBr(Test, BB_rend, BB_rjmp);

		Bld->SetInsertPoint(BB_rjmp);
		Ul = Bld->CreateCall(UlGet, CtlPtr);
		Bld->CreateBr(BB_rend);

		// common end
		Bld->SetInsertPoint(BB_rend);
		YindxAdd = Bld->CreatePHI(Type::Int64Ty, "yindx_add");
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

	Delta = ConstantInt::get(Type::Int64Ty, delta_size);

	Bld->SetInsertPoint(BB);
	Size = Bld->CreateLoad(SizePtr, "size");
	Myx0 = Bld->CreateLoad(MyxPtr, "myx0");
	Bld->CreateBr(BB_lbody);

	// Body
	Bld->SetInsertPoint(BB_lbody);
	Cnt = Bld->CreatePHI(Type::Int8Ty, "cnt");
	Myx = Bld->CreatePHI(Myx0->getType(), "myx");
	doOp(Myx);

	newMyx = Bld->CreateGEP(Myx, Delta, "new_myx");

	NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");
	Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
	Bld->CreateCondBr(Test, BB_lexit, BB_lbody);

	Cnt->addIncoming(Zero8, BB);
	Cnt->addIncoming(NextCnt, BB_lbody);

	Myx->addIncoming(Myx0, BB);
	Myx->addIncoming(newMyx, BB_lbody);

	// Exit
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

	Delta = ConstantInt::get(Type::Int64Ty, delta_size);

	Bld->SetInsertPoint(BB);
	Size = Bld->CreateLoad(SizePtr, "size");
	Yindx0 = Bld->CreateLoad(YindxPtr);
	Bld->CreateBr(BB_lbody);

	// Body
	Bld->SetInsertPoint(BB_lbody);
	Cnt = Bld->CreatePHI(Type::Int8Ty, "cnt");
	Yindx = Bld->CreatePHI(Type::Int64Ty, "yindx");

	doOp(NULL, Yindx);

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

	D = ConstantInt::get(Type::Int64Ty, delta_size);
	minusD = ConstantInt::get(Type::Int64Ty, -delta_size);

	Bld->SetInsertPoint(BB);
	Size = Bld->CreateLoad(SizePtr, "size");
	Myx0 = Bld->CreateLoad(MyxPtr, "myx0");
	Yindx0 = Bld->CreateLoad(YindxPtr);
	Bld->CreateBr(BB_lbody);

	// Body
	Bld->SetInsertPoint(BB_lbody);
	Cnt = Bld->CreatePHI(Type::Int8Ty, "cnt");
	Yindx = Bld->CreatePHI(Type::Int64Ty, "yindx");
	Myx = Bld->CreatePHI(Myx0->getType(), "myx");

	doOp(Myx, Yindx);

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

void CsxJit::BlockRowCase(BasicBlock *BB,
//                           BasicBlock *BB_lbody,
//                           BasicBlock *BB_lexit,
                          BasicBlock *BB_exit,
                          int r, int c)
{
    Value **Myx = new Value*[c+1];
    Value **Yindx = new Value*[r+1];

    Bld->SetInsertPoint(BB);

    // Elements in block-row types are stored column-wise
    Myx[0] = Bld->CreateLoad(MyxPtr);
    for (int i = 0; i < c; i++) {
        Yindx[0] = Bld->CreateLoad(YindxPtr);
        for (int j = 0; j < r; j++) {
            doOp(Myx[i], Yindx[j]);
            Yindx[j+1] = Bld->CreateAdd(Yindx[j], One64);
        }

        Myx[i+1] = Bld->CreateGEP(Myx[i], One64);
    }

//    Bld->CreateStore(Myx[c-1], MyxPtr);
    Bld->CreateBr(BB_exit);
}

void CsxJit::BlockColCase(BasicBlock *BB,
//                           BasicBlock *BB_lbody,
//                           BasicBlock *BB_lexit,
                          BasicBlock *BB_exit,
                          int r, int c)
{
    Value **Myx = new Value*[c+1];
    Value **Yindx = new Value*[r+1];

    Bld->SetInsertPoint(BB);

    // Elements in block-col types are stored row-wise
    Yindx[0] = Bld->CreateLoad(YindxPtr);
    for (int i = 0; i < r; i++) {
        Myx[0] = Bld->CreateLoad(MyxPtr);
        for (int j = 0; j < c; j++) {
            doOp(Myx[j], Yindx[i]);
            Myx[j+1] = Bld->CreateGEP(Myx[j], One64);
        }
        Yindx[i+1] = Bld->CreateAdd(Yindx[i], One64);
    }

//    Bld->CreateStore(Myx[c-1], MyxPtr);
    Bld->CreateBr(BB_exit);
}

void CsxJit::doPrint(Value *Myx, Value *Yindx)
{
	Value *X, *Xindx;
	if (Myx == NULL)
		Myx = Bld->CreateLoad(MyxPtr);

	Myx = Bld->CreatePtrToInt(Myx, Type::Int64Ty, "myx_int");
	X = Bld->CreateLoad(Xptr);
	X = Bld->CreatePtrToInt(X, Type::Int64Ty, "x_int");
	Xindx = Bld->CreateSub(Myx, X);
	Xindx = Bld->CreateAShr(Xindx, Three64);

	if (Yindx == NULL)
		Yindx = Bld->CreateLoad(YindxPtr);

	Value *V;
	V = Bld->CreateLoad(Bld->CreateLoad(Vptr));

	Bld->CreateCall3(PrintYXV, Yindx, Xindx, V);
}

void CsxJit::doMul(Value *Myx, Value *Yindx)
{
	Value *V, *X;

	if (Myx == NULL)
		Myx = Bld->CreateLoad(MyxPtr);

	X = Bld->CreateLoad(Myx, "x");
	V = Bld->CreateLoad(Bld->CreateLoad(Vptr, "v_ptr"), "val");
	V = Bld->CreateMul(V, X, "mul");

	if (Yindx == NULL){
		Value *Yr;
		// use Yr to store the result
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

void CsxJit::doIncV()
{
	Value *V = Bld->CreateLoad(Vptr);
	Value *newV = Bld->CreateGEP(V, One64);
	Bld->CreateStore(newV, Vptr);
}

void CsxJit::doOp(Value *Myx, Value *Yindx)
{
//	doPrint(Myx, Yindx);
	doMul(Myx, Yindx);
	doIncV();
}


void CsxJit::doDeltaAddMyx(int delta_bytes)
{
	Function *F;

	switch (delta_bytes){
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
	// align ctl
	if (delta_bytes > 1){
		Align = ConstantInt::get(Type::Int32Ty,delta_bytes);
		Bld->CreateCall2(AlignF, CtlPtr, Align);
	}
	Size = Bld->CreateLoad(SizePtr, "size");
	Bld->CreateBr(BB_entry);

	// Entry
	Bld->SetInsertPoint(BB_entry);

	doOp();

	Test = Bld->CreateICmpUGT(Size, One8);
	Bld->CreateCondBr(Test, BB_body, BB_exit);

	// Body
	Bld->SetInsertPoint(BB_body);
	Cnt = Bld->CreatePHI(Type::Int8Ty, "cnt");
	doDeltaAddMyx(delta_bytes);
	NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");

	doOp();

	Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
	Bld->CreateCondBr(Test, BB_exit, BB_body);

	Cnt->addIncoming(One8, BB_entry);
	Cnt->addIncoming(NextCnt, BB_body);
}

void CsxJit::doBodyHook()
{
	BasicBlock *BB, *BB_next, *BB_default, *BB_case;
	Value *PatternMask;
	Value *v;
    uint64_t delta;
    SpmIterOrder type;

	BB = llvm_hook_newbb(M, "__body_hook", SpmvF, &BB_next);

	// get pattern for switch instruction
	Bld->SetInsertPoint(BB);
	PatternMask = ConstantInt::get(Type::Int8Ty, CTL_PATTERN_MASK);
	v = Bld->CreateLoad(FlagsPtr, "flags");
	v = Bld->CreateAnd(PatternMask, v, "pattern");

	// switch default block (call the fail function)
	BB_default = BasicBlock::Create("default", BB->getParent(), BB_next);
	Bld->SetInsertPoint(BB_default);
	Bld->CreateCall(FailF);
	Bld->CreateBr(BB_next);

	// switch instruction
	SwitchInst *Switch;
	Bld->SetInsertPoint(BB);
	//std::cerr << "Constructing switch with " << CsxMg->patterns.size() << " cases\n";
	Switch = Bld->CreateSwitch(v, BB_default, CsxMg->patterns.size());

	// Fill up switch, by iterating given patterns
	CsxManager::PatMap::iterator pat_i = CsxMg->patterns.begin();
	BasicBlock *BB_lentry, *BB_lbody, *BB_lexit;
	for ( ; pat_i !=  CsxMg->patterns.end(); ++pat_i){
		//std::cerr << "pat:" << pat_i->first << " flag:" << (int)pat_i->second.flag << "\n";

		// Alocate case + loop BBs
		BB_case = BasicBlock::Create("case", BB->getParent(), BB_default);
        
        type  = static_cast<SpmIterOrder>(pat_i->first / PID_OFFSET);
        delta = pat_i->first % PID_OFFSET;
		switch (type){
            // Deltas
        case 0:
            assert(delta ==  8 ||
                   delta == 16 ||
                   delta == 32 ||
                   delta == 64);
//        case 8: case 16: case 32: case 64:
			std::cout << "type:DELTA size:" << delta << " elements:" << pat_i->second.nr << "\n";
			BB_lentry = BasicBlock::Create("lentry", BB->getParent(), BB_default);
			BB_lbody = BasicBlock::Create("lbody", BB->getParent(), BB_default);
			DeltaCase(BB_case,
			          BB_lentry, BB_lbody,
			          BB_next,
			          delta / 8);
			break;

			// Horizontal
        case HORIZONTAL:
//        case 10000 ... 19999:
//			delta = pat_i->first - 10000;
			std::cout << "type:DRLE order:HORIZONTAL delta:" << delta << " elements:" << pat_i->second.nr << "\n";
			BB_lbody = BasicBlock::Create("lbody", BB->getParent(), BB_default);
			BB_lexit = BasicBlock::Create("lexit", BB->getParent(), BB_default);
			HorizCase(BB_case,
			          BB_lbody, BB_lexit,
			          BB_next,
			          delta);
			break;

			// Vertical
        case VERTICAL:
//        case 20000 ... 29999:
//			delta = pat_i->first - 20000;
			std::cout << "type:DRLE order:VERTICAL delta:" << delta << " elements:" << pat_i->second.nr << "\n";
			BB_lbody = BasicBlock::Create("lbody", BB->getParent(), BB_default);
			VertCase(BB_case,
			         BB_lbody,
			         BB_next,
			         delta);
			break;

			// Diagonal
        case DIAGONAL:
//         case 30000 ... 39999:
// 			delta = pat_i->first - 30000;
			std::cout << "type:DRLE order:DIAGONAL delta:" << delta << " elements:" << pat_i->second.nr << "\n";
			BB_lbody = BasicBlock::Create("lbody", BB->getParent(), BB_default);
			DiagCase(BB_case,
			         BB_lbody,
			         BB_next,
			         delta,
			         false);
			break;

			// rdiag
        case REV_DIAGONAL:
//         case 40000 ... 49999:
// 			delta = pat_i->first - 40000;
			std::cout << "type:DRLE order:REV_DIAGONAL delta:" << delta << " elements:" << pat_i->second.nr << "\n";
			BB_lbody = BasicBlock::Create("lbody", BB->getParent(), BB_default);
			DiagCase(BB_case,
			         BB_lbody,
			         BB_next,
			         delta,
			         true);
			break;

        case BLOCK_TYPE_START ... BLOCK_COL_START - 1:
            // This is a block row type
 			std::cout << "type:block_row: "
                       << type - BLOCK_TYPE_START << "x" << delta
                       << " nnz: " << pat_i->second.nr << std::endl;
// 			BB_lbody = BasicBlock::Create("lbody", BB->getParent(), BB_default);
// 			BB_lexit = BasicBlock::Create("lexit", BB->getParent(), BB_default);
			BlockRowCase(BB_case, BB_next,
                         type - BLOCK_TYPE_START, delta);
            break;
        case BLOCK_COL_START ... BLOCK_TYPE_END:
            // This is a block col type
			std::cout << "type:block_col: "
                      << delta << "x" << type - BLOCK_COL_START
                      << " nnz: " << pat_i->second.nr << std::endl;
// 			BB_lbody = BasicBlock::Create("lbody", BB->getParent(), BB_default);
// 			BB_lexit = BasicBlock::Create("lexit", BB->getParent(), BB_default);
			BlockColCase(BB_case, BB_next,
                         delta, type - BLOCK_COL_START);
            break;
        default:
			assert(false);
		}

		Switch->addCase(
			ConstantInt::get(Type::Int8Ty, pat_i->second.flag),
			BB_case
		);
	}
}

void CsxJit::doHooks()
{
	doNewRowHook();
	doBodyHook();
}

void *CsxJit::doJit()
{
	ExecutionEngine *JIT;

	verifyModule(*M, AbortProcessAction, 0);
	//ModuleToFile(M, "M.llvm.bc");
	//doOptimize(M);
	//M->dump();
	//std::cerr << "Generating Function\n";
	JIT = SingleModule::getJIT(M);
	return JIT->getPointerToFunction(SpmvF);
}

void doEncode(SPM *Spm)
{
	DRLE_Manager *DrleMg;
	SpmIterOrder type;

	// 255-1 is because we need drle with <= 255-1 size, so that
	// patterns with jmps have 255 elements
	DrleMg = new DRLE_Manager(Spm, 4, 255-1);
	DrleMg->genAllStats();
	DrleMg->outStats(std::cerr);

	type = DrleMg->chooseType();
	if (type == NONE)
		return;
	std::cerr << "Encode for " << SpmTypesNames[type] << std::endl;
	Spm->Transform(type);
	DrleMg->Encode();
	Spm->Transform(HORIZONTAL);
	//Spm->Print(std::cerr);
}
