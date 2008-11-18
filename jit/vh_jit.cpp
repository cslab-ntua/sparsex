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
	char spmv_module[] = "crsvh_spmv_template.llvm.bc";

	HuffmanJit *hj = new HuffmanJit(tree, bits, type, HUFF_LLVM_TEMPLATE);
	hj->CreateSymDecoder();
	hj->CreateOutDecoder();
	LinkFileToModule(hj->M, spmv_module);
	//hj->M->dump();
	doOptimize(hj->M);

	return (void *)hj;
}

void *vhjit_get_spmvfn(void *hj_wrap, int ci_bits)
{
	HuffmanJit *hj = (HuffmanJit *)hj_wrap;

	char fnname[BUFF_SIZE];
	snprintf(fnname, BUFF_SIZE,
	         "spm_crs%d_vhjit_%s_mul_template",
	         ci_bits, type_name(hj->getOutType()) );

	Function *SpmvFn = hj->M->getFunction(fnname);
	assert(SpmvFn);
	return (void *)hj->JIT->getPointerToFunction(SpmvFn);
}

void vhjit_shut(void *hj_wrap)
{
	//HuffmanJit *hj = (HuffmanJit *)hj_wrap;
}

}
