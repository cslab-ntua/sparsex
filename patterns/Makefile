.PHONY: all clean

all: spmv draw

## Configuration
dynarray_dir = ../../dynarray


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
llvm_deps       = core analysis executionengine jit native bitreader ipo linker bitwriter
LLVM_CFLAGS     = $(shell llvm-config --cflags)
LLVM_CXXFLAGS   = $(shell llvm-config --cxxflags)
LLVM_LDFLAGS    = $(shell llvm-config --ldflags --libs $(llvm_deps))

dynarray_dep    = $(dynarray_dir)/dynarray.o
INC            += -I$(dynarray_dir)

CXX_COMPILE     = $(GXX) $(CXXFLAGS) $(INC) $(DEFS)

mmf.o: mmf.cc spm.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

spm.o: spm.cc spm.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

#draw.o: draw.cc spm.h draw.h
#	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

draw: draw.cc spm.h draw.h spm.o mmf.o
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_FLAGS) $(CAIROMM_LIBS) $(dynarray_dep) spm.o mmf.o $<  -o $@

csx.o: csx.cc csx.h spm.h delta.h
	$(CXX_COMPILE) -I$(dynarray_dir) $< -c -o $@

drle.o: drle.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

test.o: test.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

#drle_shell: mmf.o spm.o drle.o draw.o drle_shell.cc
#	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_LIBS) $^ -o $@

drle_test: mmf.o spm.o drle.o test.o drle_test.cc
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $^ -o $@

csx_llvm_tmpl.llvm.bc: csx_llvm_tmpl.c ctl_ll.h
	llvm-gcc -I../ $(INC) -Wall -emit-llvm -c csx_llvm_tmpl.c -o $@

llvm_jit_help.o: llvm_jit_help.cc llvm_jit_help.h
	$(GXX) $(CXXFLAGS) $(LLVM_CXXFLAGS) -c $<

jit.o: jit.cc jit.h spm.h drle.h llvm_jit_help.h
	$(CXX_COMPILE)  $(LLVM_CXXFLAGS) -c $<

spmv.o: spmv.cc
	$(CXX_COMPILE) $(LLVM_CXXFLAGS) -c $<

spmv: spmv.o jit.o llvm_jit_help.o mmf.o spm.o csx.o drle.o ../libspmv.o
	$(GXX) -ggdb -rdynamic -lnuma -Xlinker --allow-multiple-definition $^ $(LLVM_LDFLAGS) -o $@

clean:
	rm -f *.o drle_shell drle_test jit
