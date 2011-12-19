/*
 * main.cc -- Main program for invoking CSX.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "spmv.h"

int main(int argc, char **argv)
{   
    spm_mt_t *spm_mt;
    
    if (argc < 2){
        std::cerr << "Usage: " << argv[0] << " <mmf_file> ... \n";
        exit(1);
    }

    // Initialize the CSX JIT execution engine
    CsxExecutionEngine &engine = CsxJitInit();
    for (int i = 1; i < argc; i++) {
        std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
        spm_mt = GetSpmMt(argv[i], engine);
        CheckLoop(spm_mt, argv[i]);
        std::cerr.flush();
        BenchLoop(spm_mt, argv[i]);
        std::cout << "=== END BENCHMARK ===" << std::endl;
        PutSpmMt(spm_mt);
    }
    
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
