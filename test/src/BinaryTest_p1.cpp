/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/*
 * \file BinaryTest_p1.cpp
 * \brief Simple program for testing matrix caching feature
 *
 * \author Athena Elafrou
 * \date 2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/CsxBuild.hpp>
#include <sparsex/internals/CsxGetSet.hpp>
#include <sparsex/internals/CsxSaveRestore.hpp>
#include <sparsex/internals/CsxUtil.hpp>
#include <sparsex/internals/Mmf.hpp>
#include <sparsex/internals/SparseMatrix.hpp>
#include <sparsex/internals/Timer.hpp>
#include <sparsex/internals/Utility.hpp>
#include <sparsex/internals/logger/LoggerUtil.hpp>
#include <iostream>

using namespace std;
using namespace sparsex::csx;
using namespace sparsex::io;
using namespace sparsex::utilities;

static const char *program_name;

int main(int argc, char **argv)
{
    EnableInfoConsole();
    RuntimeConfiguration &config = RuntimeConfiguration::GetInstance();
    config.LoadFromEnv();
    config.CheckProperties<int, double>();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    rt_context.SetRuntimeContext(config);   

    program_name = argv[0];
    char c;
    while ((c = getopt(argc, argv, "bsh")) != -1) {
        switch (c) {
        case 'h':
            // PrintUsage(std::cerr);
            exit(0);
        default:
            // PrintUsage(std::cerr);
            exit(1);
        }
    }
    
    int remargc = argc - optind; // remaining arguments
    if (remargc < 1) {
        // PrintUsage(std::cerr);
        exit(1);
    }
    argv = &argv[optind];

    sparsex::timing::Timer timer;
    SparseMatrix<MMF<int, double> > matrix(argv[0]);
    spm_mt_t *spm_mt = 0;
    timer.Start();
    spm_mt = matrix.CreateCsx();
    timer.Pause();
    double pt_time = timer.ElapsedTime();

    timer.Clear();
    timer.Start();
    matrix.Save("csx_file");
    timer.Pause();
    double dump_time = timer.ElapsedTime();

    cout << "PART1 m: " << basename(argv[0]) 
         << " pt: " << pt_time 
         << " dump: " << dump_time 
         << endl;

    PutSpmMt<int, double>(spm_mt);

    return 0;
}
