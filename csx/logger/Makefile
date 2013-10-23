.PHONY: all clean

all: logger_test

common_dir = ../

## Configuration
CXX			?= g++
CXXFLAGS	?= -g -Wall -O3 -Wdisabled-optimization -fPIC
COMPILE		= $(CXX) $(CXXFLAGS)

%.o: %.cc %.hpp
	$(COMPILE) -c $<

logger_test: logger_test.o Logger.o OutputPolicy.o
	g++ -Wl,--allow-multiple-definition -Wl,--rpath=$(shell pwd) -o $@ $^

clean:
	rm -rf *.o logger_test
