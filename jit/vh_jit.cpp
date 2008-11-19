#if 0
#include <iostream>

#include "llvm/DerivedTypes.h"
#include "llvm/Bitcode/ReaderWriter.h"
#include "llvm/Linker.h"
//using namespace llvm;
#endif

extern "C" {
	#include "huffman.h"
	#include <stdio.h>
}

#include "huff_llvm_jit.h"
#include "llvm_jit_help.h"

extern "C" {

#include "vh_jit.h"

static inline const char *type_name(hj_out_t t)
{
	switch (t){
		case HJ_OUT_DOUBLE:
		return "double";

		case HJ_OUT_FLOAT:
		return "float";

		default:
		assert(0);
	}
}

#define BUFF_SIZE 128

void *vhjit_init(huff_node_t *tree, unsigned long bits, hj_out_t type)
{
	HuffmanJit *hj = new HuffmanJit(tree, bits, type, HUFF_LLVM_TEMPLATE);
	hj->CreateSymDecoder();
	hj->CreateOutDecoder();

	// link template only once
	static int vjinits = 0;
	if (vjinits++ == 0){
		const char spmv_module[] = "crsvh_spmv_template.llvm.bc";
		LinkFileToModule(hj->M, spmv_module);
	}

	//hj->M->dump();
	return (void *)hj;
}

#define TNAME "spm_crs%d_vhjit_%s_mul_template"
#define FNAME "spm_crs%d_vhjit_%s_mul_%d"
#define HNAME "__huff_decode_hook"
void *vhjit_get_spmvfn(void *hj_wrap, int ci_bits)
{
	HuffmanJit *hj = (HuffmanJit *)hj_wrap;
	const char *otype = type_name(hj->getOutType());

	// spmv template name
	char tname[BUFF_SIZE];
	snprintf(tname, BUFF_SIZE, TNAME, ci_bits, otype);

	// spmv function name
	char fname[BUFF_SIZE];
	static int cnt=0;
	snprintf(fname, BUFF_SIZE, FNAME, ci_bits, otype, cnt++);

	Module *M = hj->M;
	Function *SpmvTe = M->getFunction(tname); // template function
	Function *DecOut = hj->getOutDecodeFn();
	Function *DecSym = hj->getSymDecodeFn();
	Function *SpmvFn = CloneAndReplaceHook(M, SpmvTe, DecOut, HNAME, fname);

	// force inline
	InlineFntoFn(DecSym, DecOut);
	InlineFntoFn(DecOut, SpmvFn);
	doOptimize(M);

	M->dump();
	assert(SpmvFn);
	return (void *)hj->JIT->getPointerToFunction(SpmvFn);
}
#undef TNAME
#undef FNAME
#undef HNAME

void vhjit_shut(void *hj_wrap)
{
	//HuffmanJit *hj = (HuffmanJit *)hj_wrap;
}

}
