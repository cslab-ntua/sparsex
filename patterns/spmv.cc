#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <string>

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
    #include "../../prfcnt/timer.h"
}

#define BUFFER_SIZE 50000

using namespace csx;

typedef struct parameters {
	unsigned int thread_no;
	SPM *Spm;
	uint64_t wsize;
	char * buffer;
        int * xform_buf;
        double sampling_prob;
        uint64_t samples_max;
} Parameters;

void *thread_function(void *initial_data)
{
    Parameters *data = (Parameters *) initial_data;
    unsigned int thread_no = data->thread_no;
    SPM *Spm = data->Spm;
    uint64_t wsize = data->wsize;
    char *buffer = data->buffer;
    int *xform_buf = data->xform_buf;
    double sampling_prob = data->sampling_prob;
    uint64_t samples_max = data->samples_max;

    char number[4];
    DRLE_Manager *DrleMg;

    strcpy(buffer,"==> Thread: #");
    sprintf(number,"%d",thread_no);
    strcat(buffer,number);
    strcat(buffer,"\n");
    DrleMg = new DRLE_Manager(Spm, 4, 255-1, 0.1, wsize, DRLE_Manager::SPLIT_BY_NNZ, sampling_prob, samples_max);
    DrleMg->ignoreAll();
    int i=0;
    strcat(buffer,"Encoding type: ");
    while (xform_buf[i] != -1) {
        int t = xform_buf[i++];
        DrleMg->removeIgnore(static_cast<SpmIterOrder>(t));
        if (i != 1)
            strcat(buffer,", ");
        strcat(buffer,SpmTypesNames[t]);
    }
    strcat(buffer,"\n");
    DrleMg->EncodeAll(buffer);
    //DrleMg->EncodeSerial(xform_buf);
    delete DrleMg;
    return 0;
}

static spm_mt_t *getSpmMt(char *mmf_fname)
{
    unsigned int threads_nr, *threads_cpus;
    spm_mt_t *spm_mt;
    SPM *Spms;
    char ** buffer;
    int *xform_buf;
    Parameters *data;
    pthread_t *threads;
    CsxJit **Jits;

    mt_get_options(&threads_nr, &threads_cpus);

    spm_mt = (spm_mt_t *) malloc(sizeof(spm_mt_t));
    if (!spm_mt){
        perror("malloc");
        exit(1);
    }

    spm_mt->nr_threads = threads_nr;
    spm_mt->spm_threads = (spm_mt_thread_t *) malloc(sizeof(spm_mt_thread_t)*threads_nr);
    if (!spm_mt->spm_threads){
        perror("malloc");
        exit(1);
    }

    for (unsigned int i=0; i<threads_nr; i++)
        spm_mt->spm_threads[i].cpu = threads_cpus[i];
    Spms = SPM::loadMMF_mt(mmf_fname, threads_nr);

    const char *wsize_str = getenv("WINDOW_SIZE");
    uint64_t wsize;

    if (!wsize_str)
        wsize = 0;
    else
        wsize = atol(wsize_str);

    const char  *sampling_prob_str = getenv("SAMPLING_PROB");
    const char  *samples = getenv("SAMPLES");
    double      sampling_prob;
    uint64_t    samples_max;

    if (!sampling_prob_str)
        sampling_prob = 0.0;
    else
        sampling_prob = atof(sampling_prob_str);

    if (!samples)
        samples_max = std::numeric_limits<uint64_t>::max();
    else
        samples_max = atol(samples);

    xtimer_t timer;
    timer_init(&timer);
    timer_start(&timer);

    char *xform_orig = getenv("XFORM_CONF");
    xform_buf = (int *) malloc(XFORM_MAX*sizeof(int));
    if (!xform_buf){
        perror("malloc");
        exit(1);
    }

    if (xform_orig) {
        int t = atoi(strtok(xform_orig, ","));
	int j=0;
        char *token;

	xform_buf[j] = t;
	j++;
        while ( (token = strtok(NULL, ",")) != NULL) {
            t = atoi(token);
    	    xform_buf[j] = t;
            j++;
        }
	while (j < XFORM_MAX) {
    	    xform_buf[j] = -1;
            j++;
        }
    }

    threads = (pthread_t *) malloc((threads_nr-1)*sizeof(pthread_t));
    if (!threads){
        perror("malloc");
        exit(1);
    }

    data = (Parameters *) malloc(threads_nr*sizeof(Parameters));
    if (!data){
        perror("malloc");
        exit(1);
    }

    buffer = (char **) malloc(threads_nr*sizeof(char *));
    if (!buffer){
        perror("malloc");
        exit(1);
    }

    Jits = (CsxJit **) malloc(threads_nr*sizeof(CsxJit *));
    if (!Jits){
        perror("malloc");
        exit(1);
    }

    for (unsigned int i=0; i < threads_nr; i++) {
        buffer[i] = (char *) malloc(BUFFER_SIZE*sizeof(char));
        if (!buffer[i]){
            perror("malloc");
            exit(1);
        }
        data[i].Spm = &Spms[i];
        data[i].wsize = wsize;
        data[i].thread_no = i;
        data[i].buffer = buffer[i];
        data[i].xform_buf = xform_buf;
        data[i].sampling_prob = sampling_prob;
        data[i].samples_max = samples_max;
    }

    for (unsigned int i=1; i<threads_nr; i++)
        pthread_create(&threads[i-1],NULL,thread_function,(void *) &data[i]);
    thread_function((void *) &data[0]);
    for (unsigned int i=0; i < threads_nr; i++){
        if (i > 0)
            pthread_join(threads[i-1],NULL);
        CsxManager *CsxMg = new CsxManager(&Spms[i]);
	spm_mt->spm_threads[i].spm = CsxMg->mkCsx();
        Jits[i] = new CsxJit(CsxMg, i);
        Jits[i]->doHooks(buffer[i]);
        delete CsxMg;
        std::cout << buffer[i];
        free(buffer[i]);
    }
    doOptimize(Jits[0]->M);
    for (unsigned int i=0; i < threads_nr; i++){
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->doJit();
        delete Jits[i];
    }
    free(Jits);
    free(buffer);
    free(data);
    free(threads);
    delete[] Spms;
    timer_pause(&timer);
    std::cout << "Preprocessing time: "
              << timer_secs(&timer) << " sec" << std::endl;
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
        std::cout << basename(argv[i]) << ": " << std::endl;
        spm_mt = getSpmMt(argv[i]);
     	CheckLoop(spm_mt, argv[i]);
        std::cerr.flush();
        BenchLoop(spm_mt, argv[i]);
        putSpmMt(spm_mt);
    }

    return 0;
}
