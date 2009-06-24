.PHONY: all clean

all: drle_shell drle_test

GXX            ?= g++
CXXFLAGS       ?= -Wall -O3 -Wdisabled-optimization
#CXXFLAGS       += -Winline
CXXFLAGS       += -ggdb -rdynamic
DEFS            =
CXX_COMPILE     = $(GXX) $(CXXFLAGS) $(INC) $(DEFS)
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS   =  $(shell pkg-config cairomm-1.0 --libs)

mmf.o: mmf.cc spm.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

spm.o: spm.cc spm.h mmf.h
	$(CXX_COMPILE) $< -c -o $@

draw.o: draw.cc spm.h draw.h
	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

drle.o: drle.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

test.o: test.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

drle_shell: mmf.o spm.o drle.o draw.o drle_shell.cc
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $(CAIROMM_LIBS) $^ -o $@

drle_test: mmf.o spm.o drle.o test.o drle_test.cc
	$(CXX_COMPILE) -Xlinker --allow-multiple-definition $^ -o $@

clean:
	rm -f *.o drle_shell drle_test
