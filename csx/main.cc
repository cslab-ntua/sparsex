/*
 * main.cc -- Main program for invoking CSX.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "spmv.h"

static const char *program_name;

void PrintUsage(std::ostream &os)
{
    os << "Usage: " << program_name
       << " [-s] [-b] <mmf_file> ...\n";
}

void PrintHelp(std::ostream &os)
{
    os << "-s    Use CSX for symmetric matrices.\n"
       << "-b    Disable the split-blocks optimization.\n"
       << "-h    Print this help message and exit.\n";
}

int main(int argc, char **argv)
{   
    char c;
    bool split_blocks = true;
    bool symmetric = false;
    spm_mt_t *spm_mt;

    program_name = argv[0];
    while ((c = getopt(argc, argv, "bsh")) != -1) {
        switch (c) {
        case 'b':
            split_blocks = false;
            break;
        case 's':
            symmetric = true;
            break;
        case 'h':
            PrintHelp(std::cerr);
            exit(0);
        default:
            PrintUsage(std::cerr);
            exit(1);
        }
    }
    
    int remargc = argc - optind; // remaining arguments
    if (remargc < 1) {
        PrintUsage(std::cerr);
        exit(1);
    }
    argv = &argv[optind];
    
    // Initialize the CSX JIT execution engine
    CsxExecutionEngine &engine = CsxJitInit();
    for (int i = 0; i < remargc; i++) {    
        std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
        spm_mt = GetSpmMt(argv[i], engine, split_blocks, symmetric);
        CheckLoop(spm_mt, argv[i]);
        std::cerr.flush();
        BenchLoop(spm_mt, argv[i]);
        std::cout << "=== END BENCHMARK ===" << std::endl;
        PutSpmMt(spm_mt);
    }
    
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
