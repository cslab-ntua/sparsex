#include <iostream>

#include "llvm_jit_help.h"

namespace csx {

} // csx namespace end

int main(int argc, char **argv)
{
	Module *M;
	AnnotationMap *annotations;

	M = ModuleFromFile("ctl_llvm_tmpl.llvm.bc");
	annotations = makeAnnotationMap(M);
	dumpAnnotationMap(annotations);

	//doOptimize(M);
	//M->dump();

	return 0;
}
