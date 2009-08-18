
#include <cstdlib>

#include "spm.h"
#include "mmf.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"
#include "llvm_jit_help.h"

extern "C" {
	#include "../mt_lib.h"
	#include "../spmv_method.h"
	#include "../spm_crs.h"
	#include "../spm_mt.h"
	#include "../spmv_loops_mt.h"
}

using namespace csx;

static spm_mt_t *getSpmMt(char *mmf_fname)
{
	unsigned int threads_nr, *threads_cpus;
	spm_mt_t *spm_mt;
	SPM *Spms, *Spm;
	DRLE_Manager *DrleMg;
	CsxManager *CsxMg;
	CsxJit **Jits;

	mt_get_options(&threads_nr, &threads_cpus);
	std::cout << "Using " << threads_nr << " threads\n";

	spm_mt = (spm_mt_t *)malloc(sizeof(spm_mt_t));
	if (!spm_mt){
		perror("malloc");
		exit(1);
	}

	spm_mt->nr_threads = threads_nr;
	spm_mt->spm_threads = (spm_mt_thread_t *)malloc(sizeof(spm_mt_thread_t)*threads_nr);
	if (!spm_mt->spm_threads){
		perror("malloc");
		exit(1);
	}

	Spms = SPM::loadMMF_mt(mmf_fname, threads_nr);
	Jits = new CsxJit *[threads_nr];

	for (unsigned int i=0; i < threads_nr; i++){
		spm_mt_thread_t *spm_mt_thread;

		std::cout << "==> Thread: #" << i << "\n";

		Spm = Spms + i;
		spm_mt_thread = spm_mt->spm_threads + i;

		DrleMg = new DRLE_Manager(Spm, 4, 255-1);
		DrleMg->EncodeAll();

		CsxMg = new CsxManager(Spm);
		spm_mt_thread->spm = CsxMg->mkCsx();

		Jits[i] = new CsxJit(CsxMg, i);
		Jits[i]->doHooks();
		spm_mt_thread->cpu = threads_cpus[i];

		delete CsxMg;
		delete DrleMg;
	}

	doOptimize(Jits[0]->M);

	for (unsigned int i=0; i < threads_nr; i++){
		spm_mt->spm_threads[i].spmv_fn = Jits[i]->doJit();
		delete Jits[i];
	}

	free(threads_cpus);
	delete[] Jits;
	delete[] Spms;

	return spm_mt;
}

static void putSpmMt(spm_mt_t *spm_mt)
{
	free(spm_mt->spm_threads);
	free(spm_mt);
}

static void CheckLoop(spm_mt_t *spm_mt, char *mmf_name)
{
	void *crs;
	uint64_t nrows, ncols, nnz;

	crs = spm_crs32_double_init_mmf(mmf_name, &nrows, &ncols, &nnz);
	std::cout << "Checking ... ";
	spmv_double_check_mt_loop(crs, spm_mt,
	                          spm_crs32_double_multiply, 1,
	                          nrows, ncols,
	                          NULL);
	spm_crs32_double_destroy(crs);
	std::cout << "Check Passed\n";
}

static void BenchLoop(spm_mt_t *spm_mt, char *mmf_name)
{
	uint64_t nrows, ncols, nnz;
	double secs, flops;
	long loops_nr = 128;

	getMmfHeader(mmf_name, nrows, ncols, nnz);
	secs = spmv_double_bench_mt_loop(spm_mt, loops_nr, nrows, ncols, NULL);
	flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
	printf("m:%s f:%s s:%lu t:%lf r:%lf\n",
	        "csx", basename(mmf_name), 0UL, secs, flops);
}

int main(int argc, char **argv)
{
	spm_mt_t *spm_mt;

	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " <mmf_file>\n";
		exit(1);
	}

	spm_mt = getSpmMt(argv[1]);
	CheckLoop(spm_mt, argv[1]);
	BenchLoop(spm_mt, argv[1]);
	putSpmMt(spm_mt);

	return 0;
}
