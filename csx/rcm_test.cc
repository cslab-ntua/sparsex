#include "rcm.h"
#include "runtime.h"
#include "spmv.h"

int main(int argc, char **argv)
{
    SPM *spms = NULL;
    spm_mt_t *spm_mt;

    bool zero_based(true);
    size_t nr_rows = 10;
    size_t nr_cols = 10;
    // sym3.mm
    /*uint64_t rowptr[] = {1,4,5,8,10,12,13,16};
    uint64_t colind[] = {1,3,7,2,1,3,4,3,4,5,7,6,1,5,7};
    uint64_t rowptr[] = {0,3,4,7,9,11,12,15};
    uint64_t colind[] = {0,2,6,1,0,2,3,2,3,4,6,5,0,4,6};
    double values[] = {1,2,3,4,2,5,6,6,7,8,9,10,3,9,11};*/
    // sym2.mm
    uint64_t rowptr[] = {0,3,8,12,17,21,25,30,32,34,36};
    uint64_t colind[] = {0,3,5,1,2,4,6,9,1,2,3,4,0,2,3,5,8,1,2,4,6,0,3,5,6,1,4,5,6,7,6,7,3,8,1,9};
    double values[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1};
    // demopatt.mtx.sorted
    /*uint64_t rowptr[] = {1,6,7,11,16,19,23,25,30,34,39};
    uint64_t colind[] = {1,2,3,4,9,8,1,2,7,10,1,2,4,6,10,1,2,10,1,2,6,10,3,4,3,4,5,6,8,3,4,5,6,3,4,5,6,10};
    double values[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};*/

    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    CsxContext csx_context;
    Configuration config;
    config = ConfigFromEnv(config);
    rt_context.SetRuntimeContext(config);
    csx_context.SetCsxContext(config);

    //spms = LoadCSR_RCM<uint64_t, double>(rowptr, colind, values, nr_rows, nr_cols, zero_based, rt_context);
    spms = LoadMMF_RCM(argv[1], rt_context);
   
    std::cout << "=== BEGIN BENCHMARK ===" << std::endl;
    spm_mt = BuildCsx(spms, NULL, rt_context, csx_context);
    std::cerr.flush();
    BenchLoop(spm_mt, argv[1]);
    //double imbalance = CalcImbalance(spm_mt);
    //std::cout << "Load imbalance: " << 100*imbalance << "%\n";
    std::cout << "=== END BENCHMARK ===" << std::endl;
    PutSpmMt(spm_mt);
    
    return 0;
}
