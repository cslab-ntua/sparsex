.PHONY: all clean

all: drle_stats

GXX            ?= g++
CXXFLAGS       ?= -Wall -O3 -Wdisabled-optimization
#CXXFLAGS       += -Winline
CXXFLAGS       += -ggdb -rdynamic
DEFS            =
CXX_COMPILE     = $(GXX) $(CXXFLAGS) $(INC) $(DEFS)
CAIROMM_FLAGS   = $(shell pkg-config cairomm-1.0 --cflags)
CAIROMM_LIBS   =  $(shell pkg-config cairomm-1.0 --libs)

spm.o: spm.cc spm.h
	$(CXX_COMPILE) $< -c -o $@

draw.o: draw.cc spm.h draw.h
	$(CXX_COMPILE) $(CAIROMM_FLAGS) $< -c -o $@

drle.o: drle.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

drle_stats.o: drle_stats.cc drle.h spm.h
	$(CXX_COMPILE) $< -c -o $@

drle_stats: spm.o drle.o drle_stats.o draw.o
	$(CXX_COMPILE) $(CAIROMM_LIBS) $^ -o $@

clean:
	rm -f *.o
