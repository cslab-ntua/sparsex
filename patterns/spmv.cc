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
    #include "../../prfcnt/tsc.h"
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
	//std::cout << "Using " << threads_nr << " threads\n";
	spm_mt = (spm_mt_t *)malloc(sizeof(spm_mt_t));
	if (!spm_mt){
		perror("malloc");
		exit(1);
	}
	spm_mt->nr_threads = threads_nr;
	spm_mt->spm_threads =
        (spm_mt_thread_t *) malloc(sizeof(spm_mt_thread_t)*threads_nr);
	if (!spm_mt->spm_threads){
		perror("malloc");
		exit(1);
	}
	Spms = SPM::loadMMF_mt(mmf_fname, threads_nr);
	Jits = new CsxJit *[threads_nr];
    const char  *wsize_str = getenv("WINDOW_SIZE");
    uint64_t    wsize;
    if (!wsize_str)
        wsize = 0;
    else
        wsize = atol(wsize_str);

    const char  *sampling_prob_str = getenv("SAMPLING_PROB");
    double      sampling_prob;
    if (!sampling_prob_str)
        sampling_prob = 0.0;
    else
        sampling_prob = atof(sampling_prob_str);

    tsc_t timer;
    tsc_init(&timer);
    tsc_start(&timer);
	for (unsigned int i=0; i < threads_nr; i++){
		spm_mt_thread_t *spm_mt_thread;

		std::cout << "==> Thread: #" << i << "\n";
		Spm = Spms + i;
		spm_mt_thread = spm_mt->spm_threads + i;
		DrleMg = new DRLE_Manager(Spm, 4, 255-1, 0.1,
                                  wsize, DRLE_Manager::SPLIT_BY_NNZ,
                                  sampling_prob);
		DrleMg->ignoreAll();

		// find transformations to apply
		const char *xform_orig = getenv("XFORM_CONF");
		if (xform_orig) {
			// copy xform_orig, because the second time around wont work
			int len = strlen(xform_orig) + 1;
			char xform_string[len];
			strncpy(xform_string, xform_orig, len);
			//std::cout << xform_string << std::endl;

			int t = atoi(strtok(xform_string, ","));
			DrleMg->removeIgnore(static_cast<SpmIterOrder>(t));
			std::cout << "Encoding type: " << SpmTypesNames[t];
			char *token;
			while ( (token = strtok(NULL, ",")) != NULL) {
				t = atoi(token);
				DrleMg->removeIgnore(static_cast<SpmIterOrder>(t));
				std::cout << ", " << SpmTypesNames[t] << std::flush;
			}
			std::cout << std::endl;
		}

		//DrleMg->EncodeSerial();
		//DrleMg->MakeEncodeTree();
		DrleMg->EncodeAll();

//        Spm->PrintElems(std::cout);
//        Spm->PrintStats(std::cout);

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

    tsc_pause(&timer);
    tsc_report(&timer);
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
	std::cout << "Checking ... " << std::flush;
	spmv_double_check_mt_loop(crs, spm_mt,
	                          spm_crs32_double_multiply, 1,
	                          nrows, ncols,
	                          NULL);
	spm_crs32_double_destroy(crs);
	std::cout << "Check Passed" << std::endl << std::flush;
}

static unsigned long CsxSize(spm_mt_t *spm_mt)
{

	unsigned long ret;

	ret = 0;
	for (unsigned int i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *t = spm_mt->spm_threads + i;
		csx_double_t *csx = (csx_double_t *)t->spm;

		ret += csx->nnz*sizeof(double);
		ret += csx->ctl_size;
	}

	return ret;
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
	        "csx", basename(mmf_name), CsxSize(spm_mt), secs, flops);
}

int main(int argc, char **argv)
{
	spm_mt_t *spm_mt;

	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " <mmf_file> ... \n";
		exit(1);
	}

    for (int i = 1; i < argc; i++) {
        std::cout << basename(argv[i]) << ": " << std::endl << std::flush;
        spm_mt = getSpmMt(argv[i]);
        CheckLoop(spm_mt, argv[i]);
        std::cerr.flush();
        BenchLoop(spm_mt, argv[i]);
        putSpmMt(spm_mt);
    }

	return 0;
}
