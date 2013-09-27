.PHONY: all clean

all: output_policy.o logger.o

common_dir = ../

## Configuration
CXX			?= g++
CXXFLAGS	?= -g -Wall -O3 -Wdisabled-optimization -fPIC
COMPILE		= $(CXX) $(CXXFLAGS)

%.o: %.cc %.hpp
	$(COMPILE) -c $<

clean:
	rm -rf *.o
