#include <iostream>

#include "spm.h"
#include "ctl.h"
#include "llvm_jit_help.h"

namespace csx {

} // csx namespace end

using namespace csx;

void mkCtlJit(CtlManager *CtlMg)
{
	Module *M;
	BasicBlock *nr_BB, *nr_BB_next;
	AnnotationMap *annotations;
	IRBuilder<> *Bld;


	M = ModuleFromFile("ctl_llvm_tmpl.llvm.bc");
	Bld = new IRBuilder<>();
	annotations = makeAnnotationMap(M);
	dumpAnnotationMap(annotations);
	//doOptimize(M);
	//M->dump();

	//BB = llvm_hook_newbb(M, "__new_row_hook", &BB_next);
	//Bld->SetInsertPoint(BB);
}

int main(int argc, char **argv)
{
	SpmIdx *Spm;
	CtlManager *CtlMg;
	uint8_t *ctl;

	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " <mmf_file>\n";
		exit(1);
	}

	Spm = loadMMF_mt(argv[1], 1);
	CtlMg = new CtlManager(Spm);
	ctl = CtlMg->mkCtl();

	mkCtlJit(CtlMg);

	//delete Spm;
	delete CtlMg;

	return 0;
}
