#include "CsxBench.hpp"
#include "LoggerUtil.hpp"
#include "SparseMatrix.hpp"
#include "Types.hpp"

#include <cfloat>

static const char *program_name;

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

void PrintUsage(std::ostream &os)
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

    AlwaysUseConsole();
    UseConsole(Debug);
    //DisableInfo();

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

    SparseMatrix<MMF<spx_index_t, spx_value_t> > matrix(argv[0]);
    // matrix.Reorder();

    for (int i = 0; i < remargc; i++) {
        std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
        std::cout << "Creating CSX...\n";
        spm_mt = matrix.CreateCsx();
        CheckLoop<spx_index_t, spx_value_t>(spm_mt, argv[0]);
        std::cout << "Running 128 SpMV loops...\n";
        BenchLoop<spx_value_t>(spm_mt, argv[0]);
        // double imbalance = CalcImbalance(spm_mt);
        // std::cout << "Load imbalance: " << 100*imbalance << "%\n";
        // std::cout << "Dumping Csx to binary file...\n";
        // matrix.Save("csx_file");
        std::cout << "Convert to internal repr.: " << internal_time << std::endl;
        std::cout << "Convert to CSX: " << csx_time << std::endl;
        // std::cout << "Dump to binary file: " << dump_time << std::endl;
        std::cout << "==== END BENCHMARK ====" << std::endl;
        matrix.Destroy();
    }

    return 0;
}
