#!/bin/sh
g++ -O3 -ggdb -rdynamic -Wall spm_csx.cc  -o spm_csx $(pkg-config  cairomm-1.0 --libs --cflags)
