#include "SparseMatrix.h"
//#include "SaveRestore.h"
#include <cfloat>

static const char *program_name;

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

    // // demopatt.mtx.sorted
    // bool zero_based(true);
    // size_t nr_rows = 10;
    // size_t nr_cols = 10;
    // uint64_t rowptr[] = {1,6,7,11,16,19,23,25,30,34,39};
    // uint64_t colind[] = {1,2,3,4,9,8,1,2,7,10,1,2,4,6,10,1,2,10,1,2,6,10,3,4,3,4,5,6,8,3,4,5,6,3,4,5,6,10};
    // uint64_t rowptr[] = {0,5,6,10,15,18,22,24,29,33,38};
    // uint64_t colind[] = {0,1,2,3,8,7,0,1,6,9,0,1,3,5,9,0,1,9,0,1,5,9,2,3,2,3,4,5,7,2,3,4,5,2,3,4,5,9};
    // double values[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26.1,26.2,27,28,29,29.1,29.2,30,31,31.1,31.2,32};

    // // sym2.mm
    // uint64_t rowptr[] = {0,3,8,12,17,21,25,30,32,34,36};
    // uint64_t colind[] = {0,3,5,1,2,4,6,9,1,2,3,4,0,2,3,5,8,1,2,4,6,0,3,5,6,1,4,5,6,7,6,7,3,8,1,9};
    // double values[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1};

    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    CsxContext csx_context;
    Configuration config;
    config = ConfigFromEnv(config, symmetric, split_blocks);
    rt_context.SetRuntimeContext(config);
    csx_context.SetCsxContext(config);

    // SparseMatrix<CSR<uint64_t, double> > matrix(rowptr, colind, values, nr_rows,
    //                                             nr_cols, zero_based,
    //                                             rt_context.GetNrThreads());
    SparseMatrix<MMF<uint64_t, double> > matrix(argv[0]);
    // matrix.Reorder();

    for (int i = 0; i < remargc; i++) {    
        std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
        spm_mt = matrix.CreateCsx(rt_context, csx_context);
        //matrix.PrintEncoded(std::cout);
        //PrintCsx(spm_mt);
        CheckLoop(spm_mt, argv[i]);
        std::cerr.flush();
        BenchLoop(spm_mt, argv[i]);
        double imbalance = CalcImbalance(spm_mt);
        std::cout << "Load imbalance: " << 100*imbalance << "%\n";
        std::cout << "=== END BENCHMARK ===" << std::endl;
        PutSpmMt(spm_mt);
    }

    // // Get/Set testing
    // std::cout << matrix.GetValue(6,6) << std::endl; 
    // matrix.SetValue(6,6,24);
    // std::cout << matrix.GetValue(6,6) << std::endl;

    // // Binary testing
    // SaveCsx(spm_mt, "/home/athena/Desktop/libcsx/csx/dump");
    // csx_double_t *spm = NULL;
    // spm = RestoreCsx("/home/athena/Desktop/libcsx/csx/dump");
    // for (size_t i = 0; i < rt_context.GetNrThreads(); i++) {
    //     std::cout << "Partition " << i << std::endl;
    //     std::cout << "nnz: " << spm[i].nnz << " " << "cols:" << spm[i].ncols<< " "
    //               << "rows: " << spm[i].nrows << std::endl;
    //     std::cout << "Dump of values: ["; 
    //     for (size_t j = 0; j < spm[i].nnz; j++)
    //         std::cout << spm[i].values[j] << " ";
    //     std::cout << "]\n";
    // }
    
    return 0;
}
