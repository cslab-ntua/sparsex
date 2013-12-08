#include "../api/types.h"
#include "CsxSaveRestore.hpp"
#include "CsxGetSet.hpp"
#include "CsxBuild.hpp"
#include "CsxBench.hpp"
#include "CsxUtil.hpp"
#include "Mmf.hpp"
#include <cfloat>

static double CalcImbalance(void *arg)
{
    spm_mt_t *spm_mt = (spm_mt_t *) arg;
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

int main(int argc, char **argv)
{   
    spm_mt_t *spm_mt = 0;
    timing::Timer timer;

    std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
    std::cout << "Reconstructing CSX from binary file..." << std::endl;
    timer.Start();
    spm_mt = RestoreCsx<int, double>("csx_file", NULL);
    timer.Pause();
    double rebuild_time = timer.ElapsedTime();

    std::cout << "Assigning new values to all matrix entries..." << std::endl;
    MMF<index_t, value_t> mmf(argv[1]);
    MMF<index_t, value_t>::iterator iter = mmf.begin();
    MMF<index_t, value_t>::iterator iter_end = mmf.end();
    double assign_time;

    if (spm_mt->symmetric) {
        timer.Clear();
        timer.Start();
        for (;iter != iter_end; ++iter) {
            SetValueCsxSym<int, double>(spm_mt, (*iter).row, (*iter).col, (*iter).val);
        }
        timer.Pause();
        assign_time = timer.ElapsedTime();
    } else {
        timer.Clear();
        timer.Start();
        for (;iter != iter_end; ++iter) {
            SetValueCsx<int, double>(spm_mt, (*iter).row, (*iter).col, (*iter).val);
        }
        timer.Pause();
        assign_time = timer.ElapsedTime();
    } 

    CheckLoop<index_t, value_t>(spm_mt, argv[1]);
    std::cerr.flush();
    std::cout << "Running 128 SpMV loops..." << std::endl;
    BenchLoop<value_t>(spm_mt, argv[1]);
    double imbalance = CalcImbalance(spm_mt);
    std::cout << "Load imbalance: " << 100*imbalance << "%\n";

    std::cout << "=== END BENCHMARK ===" << std::endl;
    PutSpmMt<value_t>(spm_mt);
    std::cout << "Load from binary: " << rebuild_time << std::endl;
    std::cout << "Assign new values: " << assign_time << std::endl;

    return 0;
}
