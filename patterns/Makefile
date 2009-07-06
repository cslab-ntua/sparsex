.PHONY: all clean

all: drle_shell drle_test jit

GXX            ?= g++
CXXFLAGS       ?= -Wall -O2 -Wdisabled-optimization
#CXXFLAGS       += -Winline
CXXFLAGS       += -ggdb -rdynamic
DEFS            =
INC             =

## cairomm
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS    = $(shell pkg-config cairomm-1.0 --libs)

#llvm
llvm_deps       = core analysis executionengine jit native bitreader ipo linker
LLVM_CFLAGS     = $(shell llvm-config --cflags)
LLVM_CXXFLAGS   = $(shell llvm-config --cxxflags)
LLVM_LDFLAGS    = $(shell llvm-config --ldflags --libs $(llvm_deps))

dynarray_dir    = $(shell rsrc resource 'dynarray')
dynarray_dep    = $(dynarray_dir)/dynarray.o
INC            += -I$(dynarray_dir) 

CXX_COMPILE     = $(GXX) $(CXXFLAGS) $(INC) $(DEFS)

mmf.o: mmf.cc spm.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

spm.o: spm.cc spm.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

draw.o: draw.cc spm.h draw.h
	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

ctl.o: ctl.cc ctl.h spm.h delta.h
	$(CXX_COMPILE) -I$(dynarray_dir) $< -c -o $@

drle.o: drle.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

test.o: test.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

drle_shell: mmf.o spm.o drle.o draw.o drle_shell.cc
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_LIBS) $^ -o $@

drle_test: mmf.o spm.o drle.o test.o drle_test.cc
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $^ -o $@

ctl_llvm_tmpl.llvm.bc: ctl_llvm_tmpl.c ctl_ll.h
	llvm-gcc -Wall -emit-llvm -c ctl_llvm_tmpl.c -o $@

llvm_jit_help.o: llvm_jit_help.cc
	$(GXX) $(CXXFLAGS) $(LLVM_CXXFLAGS) -c $<

jit: jit.cc llvm_jit_help.o mmf.o spm.o ctl.o drle.o $(dynarray_dep)
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $^ $(LLVM_CXXFLAGS) $(LLVM_LDFLAGS) -o $@

clean:
	rm -f *.o drle_shell drle_test jit
