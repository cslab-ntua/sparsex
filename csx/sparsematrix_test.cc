#include "../C-API/mattype.h"
#include "sparse_matrix.h"
#include "sparse_matrix_impl.h"
#include "csx_matvec.h"
#include "csx_bench.h"
#include "logger.hpp"

#include <cfloat>

/*
 * Explicit instantiation declarations: prevent implicit instantiations.
 * Code that would otherwise cause an implicit instantiation has to
 * use the explicit instatiation definition provided somewhere else in the
 * program.
 */
extern template class SparseMatrix<MMF<index_t, value_t> >;
extern template class SparseMatrix<CSR<index_t, value_t> >;
extern template spm_mt_t *RestoreCsx<index_t, value_t>(const char *, index_t **);

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
    long loops = 128;

    // logging::AlwaysUseFile();

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

    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    CsxContext csx_context;
    Configuration config;
    config = ConfigFromEnv(config, symmetric, split_blocks);
    rt_context.SetRuntimeContext(config);
    csx_context.SetCsxContext(config);
    csx::Timer timer;
    double pt, t;

    SparseMatrix<MMF<index_t, value_t> > matrix(argv[0]);
    // matrix.Reorder();

    for (int i = 0; i < remargc; i++) {    
        std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
        if (symmetric) {
            std::cout << "Creating CSX-Sym...\n";
            spm_mt = matrix.CreateCsx<true>(csx_context, pt);
        } else {
            std::cout << "Creating CSX...\n";
            spm_mt = matrix.CreateCsx<false>(csx_context, pt);
        }
        CheckLoop(spm_mt, argv[0]);
        std::cout << "Running 128 SpMV loops..." << std::endl;
        BenchLoop(spm_mt, argv[0]);
        // vector_t *x = vec_create_random(matrix.GetNrCols());
        // vector_t *y = vec_create_random(matrix.GetNrRows());
        // if (symmetric) {
        //     timer.Start();
        //     for (int i = 0; i < loops; i++)
        //         spmv_sym_mt(spm_mt, x, 2, y, 1);
        //     timer.Pause();
        // } else {
        //     timer.Start();
        //     for (int i = 0; i < loops; i++)
        //         spmv_mt(spm_mt, x, 2, y, 1);
        //     timer.Pause();
        // }
        // t = timer.ElapsedTime();
        // std::cout << "DONE" << std::endl;
        // vec_print(y);
        // vec_destroy(x);
        // vec_destroy(y);
        // std::cout << "Dumping Csx to binary file...";
        // matrix.Save("csx_file");
        // double flops = (double)(loops*matrix.GetNrNonzeros()*2)/
        //     ((double)1000*1000*t);
        // if (symmetric)
        //     printf("m:%s ms:%lu s:%lu pt:%lf t:%lf r:%lf\n", "csx-sym",
        //            MapSize(spm_mt), CsxSymSize<double>(spm_mt), pt, t, flops);
        // else
        //     printf("m:%s s:%lu pt:%lf t:%lf r:%lf\n", "csx",
        //            CsxSize<double>(spm_mt), pt, t, flops);
        double imbalance = CalcImbalance(spm_mt);
        std::cout << "Load imbalance: " << 100*imbalance << "%\n";
        matrix.Destroy();
        std::cout << "==== END BENCHMARK ====" << std::endl;
    }

    return 0;
}
