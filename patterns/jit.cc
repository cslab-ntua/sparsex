#include <iostream>
#include <cassert>

#include "llvm/Analysis/Verifier.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ModuleProvider.h"


#include "spm.h"
#include "ctl.h"
#include "llvm_jit_help.h"


namespace csx {

} // csx namespace end

using namespace csx;
using namespace llvm;

class CtlJit {
public:
	CtlManager *CtlMg;
	Module *M;
	IRBuilder<> *Bld;
	ModuleProvider *MP;
	ExecutionEngine *JIT;
	Value *XindxPtr, *YindxPtr, *SizePtr, *FlagsPtr, *CtlPtr;

	Function  *U8GetFn;
	Function *U16GetFn;
	Function *U32GetFn;
	Function *DecodeF;
	Function *FailF;
	Function *PrintYX;
	Function *AlignF;

	Value *Zero32;
	Value *Zero64;
	Value *One8;
	Value *One32;
	Value *One64;

	Annotations annotations;

	CtlJit(CtlManager *CtlMg);
	void doNewRowHook();
	void doBodyHook();
	void doHooks();
	void *doJit();

	void DeltaCase(BasicBlock *BB,
	               BasicBlock *BB_def,
	               BasicBlock *BB_exit,
	               int delta_bytes);
};

CtlJit::CtlJit(CtlManager *_ctl_mg) : CtlMg(_ctl_mg)
{
	this->M = ModuleFromFile("ctl_llvm_tmpl.llvm.bc");
	this->Bld = new IRBuilder<>();

	this->annotations.update(M);

	this->DecodeF = M->getFunction("ctl_decode_template");
	assert(this->DecodeF);
	this->PrintYX = M->getFunction("print_yx");
	assert(this->PrintYX);
	this->FailF = M->getFunction("fail");
	assert(this->FailF);
	this->AlignF = M->getFunction("align_ptr");
	assert(this->AlignF);

	this->XindxPtr = annotations.getValue("vars::x_indx");
	this->YindxPtr = annotations.getValue("vars::y_indx");
	this->CtlPtr = annotations.getValue("vars::ctl");
	this->SizePtr = annotations.getValue("vars::size");
	this->FlagsPtr = annotations.getValue("vars::flags");

	this->Zero32 = ConstantInt::get(Type::Int32Ty, 0);
	this->Zero64 = ConstantInt::get(Type::Int64Ty, 0);
	this->One8 = ConstantInt::get(Type::Int8Ty, 1);
	this->One32 = ConstantInt::get(Type::Int32Ty, 1);
	this->One64 = ConstantInt::get(Type::Int64Ty, 1);
}

void CtlJit::doNewRowHook()
{
	BasicBlock *BB, *BB_next;
	Value *v;

	// new row
	BB = llvm_hook_newbb(M, "__new_row_hook", &BB_next);
	Bld->SetInsertPoint(BB);
	if (!CtlMg->row_jmps){
		v = Bld->CreateLoad(YindxPtr, "y_indx_load");
		v = Bld->CreateAdd(v, One64, "y_indx_inc");
		v = Bld->CreateStore(v, YindxPtr);
	} else {
		assert(false);
	}
	Bld->CreateStore(Zero64, XindxPtr);
	Bld->CreateBr(BB_next);
}

void CtlJit::DeltaCase(BasicBlock *BB, BasicBlock *BB_def, BasicBlock *BB_exit,
                       int delta_bytes)
{
	Value *Align, *Size, *Xindx, *XindxAdd, *Test, *NextCnt;
	PHINode *Cnt;
	Function *F;
	BasicBlock *BB_entry, *BB_body;
	char buff[32];

	BB_entry = BasicBlock::Create("entry", BB->getParent(), BB_def);
	BB_body = BasicBlock::Create("body", BB->getParent(), BB_def);

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
	Bld->CreateCall2(PrintYX, YindxPtr, XindxPtr);
	Test = Bld->CreateICmpUGT(Size, One8);
	Bld->CreateCondBr(Test, BB_body, BB_exit);

	// Body
	Bld->SetInsertPoint(BB_body);
	Cnt = Bld->CreatePHI(Type::Int8Ty, "cnt");
	// add delta to xIndx
	snprintf(buff, sizeof(buff), "u%d_get", delta_bytes*8);
	F = M->getFunction(buff);
	assert(F);
	Xindx = Bld->CreateLoad(XindxPtr);
	XindxAdd = Bld->CreateCall(F, CtlPtr);
	XindxAdd = Bld->CreateAdd(Xindx, XindxAdd);
	Bld->CreateStore(XindxAdd, XindxPtr);

	NextCnt = Bld->CreateAdd(Cnt, One8, "next_cnt");
	Test = Bld->CreateICmpEQ(NextCnt, Size, "cnt_test");
	Bld->CreateCall2(PrintYX, YindxPtr, XindxPtr);
	Bld->CreateCondBr(Test, BB_exit, BB_body);

	Cnt->addIncoming(One8, BB_entry);
	Cnt->addIncoming(NextCnt, BB_body);
}

void CtlJit::doBodyHook()
{
	BasicBlock *BB, *BB_next, *BB_default, *BB_case;
	Value *PatternMask;
	Value *v;

	BB = llvm_hook_newbb(M, "__body_hook", &BB_next);

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
	Switch = Bld->CreateSwitch(v, BB_default, CtlMg->patterns.size());

	// Fill up switch, by iterating given patterns
	CtlManager::PatMap::iterator pat_i = CtlMg->patterns.begin();
	for ( ; pat_i !=  CtlMg->patterns.end(); ++pat_i){

		std::cout << "pat:" << pat_i->first << " flag:" << (int)pat_i->second.flag << "\n";
		BB_case = BasicBlock::Create("case", BB->getParent(), BB_default);
		switch (pat_i->first){
			// Deltas
			case 8: case 16: case 32: case 64:
			DeltaCase(BB_case, BB_default, BB_next, pat_i->first/8);
			break;

			break;
			//horizontal
			case 10000 ... 19999:
			assert(false);
			break;
			// vertical
			case 20000 ... 29999:
			assert(false);
			break;
			// diagonal
			case 30000 ... 39999:
			assert(false);
			break;
			// rdiag
			case 40000 ... 49999:
			assert(false);
			break;
		}
		Switch->addCase(
			ConstantInt::get(Type::Int8Ty, pat_i->second.flag),
			BB_case
		);
	}
}

void CtlJit::doHooks()
{
	doNewRowHook();
	doBodyHook();
}

void *CtlJit::doJit()
{
	verifyModule(*M, AbortProcessAction, 0);
	std::cout << "Generating Function\n";
	std::string Error;
	MP = new  ExistingModuleProvider(M);
	JIT = ExecutionEngine::createJIT(MP, &Error);
	if (!JIT){
		std::cerr << "ExectionEngine::createJIT:" << Error << "\n";
		exit(1);
	}
	return JIT->getPointerToFunction(DecodeF);
}

typedef void decode_fn_t(uint8_t *ctl, unsigned long ctl_size);

int main(int argc, char **argv)
{
	SpmIdx *Spm;
	CtlManager *CtlMg;
	CtlJit *Jit;
	uint8_t *ctl;
	uint64_t ctl_size;

	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " <mmf_file>\n";
		exit(1);
	}

	Spm = loadMMF_mt(argv[1], 1);
	CtlMg = new CtlManager(Spm);
	ctl = CtlMg->mkCtl(&ctl_size);
	Jit = new CtlJit(CtlMg);

	Jit->doHooks();
	decode_fn_t *fn = (decode_fn_t *)Jit->doJit();
	fn(ctl, ctl_size);

	//delete Spm;
	delete CtlMg;

	return 0;
}
