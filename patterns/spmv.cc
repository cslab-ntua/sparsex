
#include <cstdlib>

#include "spm.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"

extern "C" {
	#include "../mt_lib.h"
	#include "../spmv_method.h"
	#include "../spm_crs.h"
	#include "../spm_mt.h"
	#include "../spmv_loops_mt.h"
}

using namespace csx;

int main(int argc, char **argv)
{
	SPM *Spms, *Spm;
	DRLE_Manager *DrleMg;
	CsxManager *CsxMg;
	CsxJit *Jit = NULL;
	csx_double_t *csx;
	spm_mt_t spm_mt;
	spm_mt_thread_t *spm_mt_thread;
	unsigned int threads_nr, *threads_cpus;

	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " <mmf_file>\n";
		exit(1);
	}

	mt_get_options(&threads_nr, &threads_cpus);
	std::cout << "Using " << threads_nr << " threads\n";

	Spms = SPM::loadMMF_mt(argv[1], threads_nr);

	spm_mt.nr_threads = threads_nr;
	spm_mt.spm_threads = (spm_mt_thread_t *)malloc(sizeof(spm_mt_thread_t)*threads_nr);
	if (!spm_mt.spm_threads){
		perror("malloc");
		exit(1);
	}

	for (unsigned int i=0; i < threads_nr; i++){
		std::cout << "************ Thread: " << i << "\n";

		Spm = Spms + i;
		spm_mt_thread = spm_mt.spm_threads + i;

		DrleMg = new DRLE_Manager(Spm, 4, 255-1);
		DrleMg->EncodeAll();
		CsxMg = new CsxManager(Spm);
		Jit = new CsxJit(CsxMg, i);

		spm_mt_thread->spm = csx = CsxMg->mkCsx();
		Jit->doHooks();
		spm_mt_thread->spmv_fn = Jit->doJit();
		spm_mt_thread->cpu = threads_cpus[i];

		delete Jit;
		delete CsxMg;
		delete DrleMg;
	}

	uint64_t nrows, ncols, nnz;
	void *crs = spm_crs32_double_init_mmf(argv[1], &nrows, &ncols, &nnz);
	spmv_double_check_mt_loop(crs, &spm_mt,
	                          spm_crs32_double_multiply, 1,
	                          nrows, ncols,
	                          NULL);
	spm_crs32_double_destroy(crs);
	std::cout << "Check Passed\n";

	free(threads_cpus);
	free(spm_mt.spm_threads);
	delete[] Spms;

	return 0;
}
