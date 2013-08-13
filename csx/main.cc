/*
 * main.cc -- Main program for invoking CSX.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include "../C-API/mattype.h"
#include "MatrixLoading.h"
#include "SparseInternal.h"
#include "CsxGetSet.h"
#include "CsxUtil.h"
#include "runtime.h"
#include "spmv.h"
#include <cstdio>
#include <cfloat>

extern template double GetValueCsx<index_t, value_t>(void *, index_t, index_t);
extern template bool SetValueCsx<index_t, value_t>(void *, index_t, index_t, value_t);
extern template void PutSpmMt<index_t, value_t>(spm_mt_t *spm_mt);

static const char *program_name;

//
// FIXME: This is a duplicate of calc_imbalance() in libspmv
// 
static double CalcImbalance(void *m)
{
    spm_mt_t *spm_mt = (spm_mt_t *) m;
    size_t i;
    
    double min_time = DBL_MAX;
    double max_time = 0.0;
    double total_time = 0.0;
    size_t worst = -1;
    for (i = 0; i < spm_mt->nr_threads; ++i) {
        spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
        double thread_time = spm->secs;
        printf("thread %zd: %f\n", i, thread_time);
        total_time += thread_time;
        if (thread_time > max_time) {
            max_time = thread_time;
            worst = i;
        }

        if (thread_time < min_time)
            min_time = thread_time;
    }

    double ideal_time = total_time / spm_mt->nr_threads;
    printf("Worst thread: %zd\n", worst);
    printf("Expected perf. improvement: %.2f %%\n",
           100*(max_time / ideal_time - 1));
    return (max_time - min_time) / min_time;
}

void PrintUsage(std::ostream &os)
{
    os << "Usage: " << program_name
       << " [-s] [-b] <mmf_file> ...\n"
       << "\t-s    Use CSX for symmetric matrices.\n"
       << "\t-b    Disable the split-blocks optimization.\n"
       << "\t-h    Print this help message and exit.\n";
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
            PrintUsage(std::cerr);
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

    // Initialization of runtime configuration
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    CsxContext csx_context;
    
    Configuration config;
    config = ConfigFromEnv(config, symmetric, split_blocks);

    rt_context.SetRuntimeContext(config);
    csx_context.SetCsxContext(config);

    // Load matrix    
    SparseInternal<index_t, value_t> *spms = NULL;
    SparsePartitionSym<index_t, value_t> *spms_sym = NULL;

    if (!csx_context.IsSymmetric()) {
        spms = LoadMMF_mt<index_t, value_t>(argv[0], rt_context.GetNrThreads());
    } else {
        spms_sym = LoadMMF_Sym_mt<index_t, value_t>(argv[0], rt_context.GetNrThreads());
    }

    for (int i = 0; i < remargc; i++) {    
        std::cout << "=== BEGIN BENCHMARK ===" << std::endl; 
        double pre_time;
        if (!symmetric)
            spm_mt = BuildCsx(spms, rt_context, csx_context, pre_time);
        else 
            spm_mt = BuildCsxSym(spms_sym, rt_context, csx_context, pre_time);
        CheckLoop(spm_mt, argv[i]);
        std::cerr.flush();
        BenchLoop(spm_mt, argv[i]);
        double imbalance = CalcImbalance(spm_mt);
        std::cout << "Load imbalance: " << 100*imbalance << "%\n";
        std::cout << "=== END BENCHMARK ===" << std::endl;

        /* Get/Set testing */
        spms = LoadMMF_mt<index_t, value_t>(argv[0], rt_context.GetNrThreads());
        SparsePartition<index_t, value_t>::PntIter iter, iter_end;      
        for (size_t i = 0; i < rt_context.GetNrThreads(); i++) {
            iter = spms->GetPartition(i)->PointsBegin();
            iter_end = spms->GetPartition(i)->PointsEnd();
            for (;iter != iter_end; ++iter) {
                index_t row = (*iter).row + spms->GetPartition(i)->GetRowStart();
                index_t col = (*iter).col;
                //std::cout << "looking for " << row << " " << col << std::endl;
                value_t value = GetValueCsx<index_t, value_t>(spm_mt, row, col);
                //std::cout << value << "=" << (*iter).val << std::endl;
                assert((*iter).val == value);
            }
        }
        /* Cleanup */
        PutSpmMt<index_t, value_t>(spm_mt);
    }
    
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
