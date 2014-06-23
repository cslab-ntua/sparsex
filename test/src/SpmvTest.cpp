#include "sparsex/internals/CsxBench.hpp"
#include "sparsex/internals/logger/LoggerUtil.hpp"
#include "sparsex/internals/SparseMatrix.hpp"
#include "sparsex/internals/Types.hpp"

#include <cfloat>

using namespace std;

static const char *program_name;

static double CalcImbalance(void *arg) __attribute__ ((unused));

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
        printf("thread %lu: %f\n", i, thread_time);
        total_time += thread_time;
        if (thread_time > max_time) {
            max_time = thread_time;
            worst = i;
        }

        if (thread_time < min_time)
            min_time = thread_time;
    }

    double ideal_time = total_time / spm_mt->nr_threads;
    printf("Worst thread: %lu\n", worst);
    printf("Expected perf. improvement: %.2f %%\n",
           100*(max_time / ideal_time - 1));
    return (max_time - min_time) / min_time;
}

void PrintUsage(ostream &os)
{
    os << "Usage: " << program_name
       << " [-s] <mmf_file> ...\n"
       << "\t-s    Use CSX for symmetric matrices.\n"
       << "\t-h    Print this help message and exit.\n";
}

int main(int argc, char **argv)
{
    char c;
    spm_mt_t *spm_mt;

#if SPX_DEBUG
    AlwaysUseConsole();
#else
    DisableInfo();
#endif

    RuntimeConfiguration &config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    config.LoadFromEnv();
    rt_context.SetRuntimeContext(config);   

    program_name = argv[0];
    while ((c = getopt(argc, argv, "bsh")) != -1) {
        switch (c) {
        case 's':
            config.SetProperty(RuntimeConfiguration::MatrixSymmetric, "true");
            break;
        case 'h':
            PrintUsage(cerr);
            exit(0);
        default:
            PrintUsage(cerr);
            exit(1);
        }
    }
    
    int remargc = argc - optind; // remaining arguments
    if (remargc < 1) {
        PrintUsage(cerr);
        exit(1);
    }
    argv = &argv[optind];

    SparseMatrix<MMF<int, double> > matrix(argv[0]);
    // matrix.Reorder();

    for (int i = 0; i < remargc; i++) {
        cout << "=== BEGIN BENCHMARK ===" << endl;
        cout << "Creating CSX...\n";
        spm_mt = matrix.CreateCsx();
        CheckLoop<unsigned int, double>(spm_mt, argv[0]);
        cout << "Running 128 SpMV loops...\n";
        BenchLoop<double>(spm_mt, argv[0]);
        // double imbalance = CalcImbalance(spm_mt);
        // cout << "Load imbalance: " << 100*imbalance << "%\n";
        // cout << "Dumping Csx to binary file...\n";
        // matrix.Save("csx_file");
        cout << "Convert to internal repr.: " << internal_time << endl;
        cout << "Convert to CSX: " << csx_time << endl;
        // cout << "Dump to binary file: " << dump_time << endl;
        cout << "==== END BENCHMARK ====" << endl;
        matrix.Destroy();
    }

    return 0;
}
