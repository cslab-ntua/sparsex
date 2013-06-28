#include "SparseInternal.h"
#include "rcm.h"
#include "runtime.h"
#include "spmv.h"

int main(int argc, char **argv)
{
    SparseInternal<uint64_t, double> *spi = NULL;
    spm_mt_t *spm_mt;

    // bool zero_based(true);
    // size_t nr_rows = 7;
    // size_t nr_cols = 7;

    // sym3.mm
    // uint64_t rowptr[] = {1,4,5,8,10,12,13,16};
    // uint64_t colind[] = {1,3,7,2,1,3,4,3,4,5,7,6,1,5,7};
    // uint64_t rowptr[] = {0,3,4,7,9,11,12,15};
    // uint64_t colind[] = {0,2,6,1,0,2,3,2,3,4,6,5,0,4,6};
    // double values[] = {1,2,3,4,2,5,6,6,7,8,9,10,3,9,11};

    // sym2.mm
    /*uint64_t rowptr[] = {0,3,8,12,17,21,25,30,32,34,36};
    uint64_t colind[] = {0,3,5,1,2,4,6,9,1,2,3,4,0,2,3,5,8,1,2,4,6,0,3,5,6,1,4,5,6,7,6,7,3,8,1,9};
    double values[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1};*/

    // demopatt.mtx.sorted
    //uint64_t rowptr[] = {1,6,7,11,16,19,23,25,30,34,39};
    //uint64_t colind[] = {1,2,3,4,9,8,1,2,7,10,1,2,4,6,10,1,2,10,1,2,6,10,3,4,3,4,5,6,8,3,4,5,6,3,4,5,6,10};
    /*uint64_t rowptr[] = {0,5,6,10,15,18,22,24,29,33,38};
    uint64_t colind[] = {0,1,2,3,8,7,0,1,6,9,0,1,3,5,9,0,1,9,0,1,5,9,2,3,2,3,4,5,7,2,3,4,5,2,3,4,5,9};*/
    //double values[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    //test4.mm
    /*uint64_t rowptr[] = {1,3,3,5,8,8,8,11,11};
    uint64_t colind[] = {1,4,2,5,2,3,6,2,5,7};
    double values[] = {1,2,3,4,8,1,6,10,9,7};*/

    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    CsxContext csx_context;
    Configuration config;
    config = ConfigFromEnv(config, false, false);
    rt_context.SetRuntimeContext(config);
    csx_context.SetCsxContext(config);

    //spi = LoadCSR_RCM_mt<uint64_t, double>(rowptr, colind, values, nr_rows,
    //                                      nr_cols, zero_based, rt_context.GetNrThreads());
    spi = LoadMMF_RCM_mt<uint64_t, double>(argv[1], rt_context.GetNrThreads());

    std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
    spm_mt = BuildCsx(spi, rt_context, csx_context);
    std::cerr.flush();
    BenchLoop(spm_mt, argv[1]);
    //double imbalance = CalcImbalance(spm_mt);
    //std::cout << "Load imbalance: " << 100*imbalance << "%\n";
    std::cout << "=== END BENCHMARK ===" << std::endl;
    PutSpmMt(spm_mt);
    
    return 0;
}
