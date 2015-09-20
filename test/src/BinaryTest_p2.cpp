/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/*
 * \file BinaryTest_p2.cpp
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
#include <sparsex/internals/Timer.hpp>
#include <sparsex/internals/Utility.hpp>
#include <iostream>

using namespace std;
using namespace sparsex::csx;
using namespace sparsex::io;
using namespace sparsex::utilities;

static const char *program_name;

int main(int argc, char **argv)
{   
    program_name = argv[0];
    spm_mt_t *spm_mt = 0;
    sparsex::timing::Timer timer;

    cout << "=== BEGIN BENCHMARK ===" << endl;
    cout << "Reconstructing CSX from binary file..." << endl;
    timer.Start();
    spm_mt = RestoreCsx<int, double>("csx_file", NULL);
    timer.Pause();
    double rebuild_time = timer.ElapsedTime();

    cout << "Assigning new values to all matrix entries..." << endl;
    MMF<int, double> mmf(argv[1]);
    int count  = mmf.GetNrNonzeros()/spm_mt->nr_threads;
    MMF<int, double>::iterator iter = mmf.begin();
    MMF<int, double>::iterator iter_end = mmf.end();
    timer.Clear();
    timer.Start();
    for (int i = 0; iter != iter_end, i < count; ++iter, i++) {
        SetValueCsx<int, double>(spm_mt, (*iter).GetRow(),
                                 (*iter).GetCol(), 0.42);
    }

    timer.Pause();
    double assign_time = timer.ElapsedTime();
    cout << "=== END BENCHMARK ===" << endl;

    PutSpmMt<int, double>(spm_mt);

    cout << "PART2 m: " << basename(argv[1]) 
         << " rebuild: " << rebuild_time 
         << " assign: " << assign_time 
         << endl;

    return 0;
}
