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

int main(int argc, char **argv)
{   
    char c;
    bool split_blocks = false;
    bool symmetric = false;
    spm_mt_t *spm_mt;

    // Initialize the CSX JIT execution engine
    CsxExecutionEngine &engine = CsxJitInit();
    
    while ((c = getopt(argc, argv, "bs")) != -1) {
        switch (c) {
        case 'b':
            split_blocks = true;
            break;

        case 's':
            symmetric = true;
            break;

        default:
            std::cerr << "Usage: " << argv[0] <<
                         " [-c] [-b] <mmf_file_1> ... \n";
            exit(1);
        }
    }
    
    int remargc = argc - optind; // remaining arguments
    
    if (remargc < 1) {
        std::cerr << "Usage: " << argv[0] << " [-c] [-b] <mmf_file_1> ... \n";
        exit(1);
    }
    argv = &argv[optind];
    
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
