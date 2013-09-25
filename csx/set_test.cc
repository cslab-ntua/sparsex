#include "../C-API/mattype.h"
#include "csx_save_restore.h"
#include "csx_get_set.h"
#include "csx_build.h"
#include "csx_bench.h"
#include "csx_util.h"
#include "mmf.h"
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
    csx::Timer timer;

    std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
    timer.Start();
    std::cout << "Reconstructing CSX/CSX-Sym from binary file..." << std::endl;
    spm_mt = RestoreCsx<int, double>("csx_file", NULL);
    std::cout << "Assigning new values to all matrix entries...";
    MMF<index_t, value_t> mmf(argv[1]);
    MMF<index_t, value_t>::iterator iter = mmf.begin();
    MMF<index_t, value_t>::iterator iter_end = mmf.end();
    if (spm_mt->symmetric) {
        for (;iter != iter_end; ++iter) {
            SetValueCsxSym<int, double>(spm_mt, (*iter).row, (*iter).col, (*iter).val);
        }
    } else {
        for (;iter != iter_end; ++iter) {
            SetValueCsx<int, double>(spm_mt, (*iter).row, (*iter).col, (*iter).val);
        }
    } 
    timer.Pause();
    double pre_time = timer.ElapsedTime();
    std::cout << "DONE" << std::endl;
    std::cout << "Total benchmark time (rebuild + assign): " << pre_time << std::endl;
    //CheckLoop(spm_mt, argv[1]);
    std::cerr.flush();
    std::cout << "Running 128 SpMV loops..." << std::endl;
    BenchLoop(spm_mt, argv[1]);
    double imbalance = CalcImbalance(spm_mt);
    std::cout << "Load imbalance: " << 100*imbalance << "%\n";
    std::cout << "=== END BENCHMARK ===" << std::endl;
    PutSpmMt<value_t>(spm_mt);

    return 0;
}
