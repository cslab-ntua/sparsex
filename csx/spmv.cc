#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <sstream>
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

#define DELTAS_MAX  64  ///> Max deltas that an encoding type may have.

using namespace csx;
using namespace std;

///> Necessary Data passed to each thread for parallel preprocessing.
typedef struct parameters {
    unsigned int    thread_no;
    SPM *           Spm;
    uint64_t        wsize;
    ostringstream   buffer;
    int             *xform_buf;
    double          sampling_prob;
    uint64_t        samples_max;
    bool            split_blocks;
    int             **deltas;
} Parameters;

int convert_string_to_string_sequence(char *str, char **str_buf,
                                      const char *start_sep,
                                      const char *end_sep)
{
    char *token = strtok(str, start_sep);
    int next = 0;
    int str_length = strcspn(token, end_sep);
    
    str_buf[next] = (char *) malloc((str_length + 1) * sizeof(char));
    strncpy(str_buf[next], token, str_length);
    str_buf[next][str_length] = 0;
    ++next;
    while ((token = strtok(NULL, start_sep)) != NULL) { 
        str_length = strcspn(token, end_sep);
        str_buf[next] = (char *) malloc((str_length - 1) * sizeof(char));
        strncpy(str_buf[next], token, str_length);
        str_buf[next][str_length] = 0;
        ++next;
    }
    
    return next;
}

void x_get_options(int **xform_buf)
{
    char *xform_orig = getenv("XFORM_CONF");
    
    *xform_buf = (int *) malloc(XFORM_MAX * sizeof(int));
    if (!(*xform_buf)) {
        perror("malloc");
        exit(1);
    }
    
    if (xform_orig && strlen(xform_orig)) {
        int next = 0;
        int t = atoi(strtok(xform_orig, ","));
        char *token;
        
        (*xform_buf)[next] = t;
        ++next;
        while ((token = strtok(NULL, ",")) != NULL) {
            t = atoi(token);
            (*xform_buf)[next] = t;
            ++next;
        }
        
        (*xform_buf)[next] = -1;
    } else {
        (*xform_buf)[0] = 0;
        (*xform_buf)[1] = -1;
    }
    
    std::cout << "Encoding type: ";
    for (unsigned int i = 0; (*xform_buf)[i] != -1; i++) {
        if (i != 0)
            std::cout << ", ";
        
        std::cout << SpmTypesNames[(*xform_buf)[i]];
    }
    
    std::cout << std::endl;
}

void encode_deltas_get_options(int ***deltas)
{
    char *encode_deltas_str = getenv("ENCODE_DELTAS");
    
    if (encode_deltas_str) {
        ///> Init matrix deltas.
        *deltas = (int **) malloc(XFORM_MAX * sizeof(int *));
        if (!*deltas) {
            perror("malloc");
            exit(1);
        }
        
        for (int i = 0; i < XFORM_MAX; i++) {
            (*deltas)[i] = (int *) malloc(DELTAS_MAX * sizeof(int));
            if (!(*deltas)[i]) {
                perror("malloc");
                exit(1);
            }
        }

        ///> Fill deltas with the appropriate data.
        char **temp = (char **) malloc(XFORM_MAX * sizeof(char *));
        int temp_size = convert_string_to_string_sequence(encode_deltas_str,
                                                          temp, "{", "}");
        
        for (int i = 0; i < temp_size; i++) {
            int j=0;
            char *token = strtok(temp[i], ",");
            
            (*deltas)[i][j] = atoi(token);
            ++j;
            while ((token = strtok(NULL, ",")) != NULL) {
                (*deltas)[i][j] = atoi(token);
                ++j;
            }
            
            (*deltas)[i][j] = -1;
        }
        
        free(temp);
        
        ///> Print deltas.
        std::cout << "Deltas to Encode: ";
        for (int i = 0; i < temp_size; i++) {
            if (i != 0)
                std::cout << "}, ";
            std::cout << "{";
            assert((*deltas)[i][0] != -1);
            std::cout << (*deltas)[i][0];
            for (int j = 1; (*deltas)[i][j] != -1; j++)
                std::cout << "," << (*deltas)[i][j];
        }
        
        std::cout << "}" << std::endl;
    }
}

uint64_t wsize_get_options()
{
    const char *wsize_str = getenv("WINDOW_SIZE");
    uint64_t wsize;
    
    if (!wsize_str) {
        wsize = 0;
        std::cout << "Window size: Not set" << std::endl;
    }
    else {
        wsize = atol(wsize_str);
        std::cout << "Window size: " << wsize << std::endl;
    }
    
    return wsize;
}

uint64_t samples_get_options()
{
    const char *samples = getenv("SAMPLES");
    uint64_t samples_max;
    
    if (!samples) {
        samples_max = std::numeric_limits<uint64_t>::max();
        std::cout << "Number of samples: Not set" << std::endl;
    }
    else {
        samples_max = atol(samples);
        std::cout << "Number of samples: " << samples_max << std::endl;
    }
    
    return samples_max;
}

double probability_get_options()
{
    const char *sampling_prob_str = getenv("SAMPLING_PROB");
    double sampling_prob;
    
    if (!sampling_prob_str) {
        sampling_prob = 0.0;
        std::cout << "Sampling prob: Not set" << std::endl;
    }
    else {
        sampling_prob = atof(sampling_prob_str);
        std::cout << "Sampling prob: " << sampling_prob << std::endl;
    }
    
    return sampling_prob;
}

bool split_get_options()
{
    const char *split_blocks_str = getenv("SPLIT_BLOCKS");
    bool split_blocks;
    
    if (!split_blocks_str)
        split_blocks = false;
    else
        split_blocks = true;
        
    return split_blocks;
}

/**
 *  Parallel Preprocessing.
 *
 *  @param initial_data parameters of matrix needed by each thread.
 */
void *thread_function(void *initial_data)
{
    Parameters *data = (Parameters *) initial_data;
    DRLE_Manager *DrleMg;

    ///> Print Thread_Id.
    data->buffer << "==> Thread: #" << data->thread_no << endl;
    /**
     *  Init DrleMg which takes statistical data and encodes and decodes the 
     *  matrix.
     */
    DrleMg = new DRLE_Manager(data->Spm, 4, 255-1, 0.1, data->wsize,
                              DRLE_Manager::SPLIT_BY_NNZ, data->sampling_prob,
                              data->samples_max);
    ///> Adjust the ignore settings properly.
    DrleMg->IgnoreAll();
    for (int i = 0; data->xform_buf[i] != -1; ++i)
        DrleMg->RemoveIgnore(static_cast<SpmIterOrder>(data->xform_buf[i]));
        
    /**
     *  If deltas choices given encode the matrix with the order given by 
     *  XFORM_CONF, else find statistical data for the types in XFORM_CONF,
     *  choose the best choise, encode it and proceed likewise until there is
     *  no satisfying encoding
     */
    if (data->deltas)
        DrleMg->EncodeSerial(data->xform_buf, data->deltas, data->split_blocks);
    else
        DrleMg->EncodeAll(data->buffer, data->split_blocks);
        
    //DrleMg->MakeEncodeTree(data->split_blocks);
    delete DrleMg;
    return 0;
}

static spm_mt_t *GetSpmMt(char *mmf_fname)
{
    unsigned int nr_threads, *threads_cpus;
    spm_mt_t *spm_mt;
    SPM *Spms;
    int *xform_buf = NULL;
    int **deltas = NULL;
    Parameters *data;
    pthread_t *threads;
    CsxJit **Jits;

    ///> Take MT_CONF.
    // TODO: In case of MT_CONF=empty a default message is printed by mt_get_options
    mt_get_options(&nr_threads, &threads_cpus);
    std::cout << "MT_CONF=";
    for (unsigned int i = 0; i < nr_threads; i++) {
        if (i != 0)
            std::cout << ",";
            
        std::cout << threads_cpus[i];
    }
    
    std::cout << std::endl;

    ///> Take XFORM_CONF.
    x_get_options(&xform_buf);

    ///> Take ENCODE_DELTAS
    encode_deltas_get_options(&deltas);

    ///> Take WINDOW_SIZE.
    uint64_t wsize = wsize_get_options();

    ///> Take SAMPLES.
    uint64_t samples_max = samples_get_options();

    ///> Take SAMPLING_PROB.
    double sampling_prob = probability_get_options();

    ///> Take SPLIT_BLOCKS
    bool split_blocks = split_get_options();

    ///> Initalization of spm_mt.
    spm_mt = (spm_mt_t *) malloc(sizeof(spm_mt_t));
    if (!spm_mt) {
        perror("malloc");
        exit(1);
    }
    
    spm_mt->nr_threads = nr_threads;
    spm_mt->spm_threads =
        (spm_mt_thread_t *) malloc(sizeof(spm_mt_thread_t) * nr_threads);
    if (!spm_mt->spm_threads){
        perror("malloc");
        exit(1);
    }
    
    for (unsigned int i = 0; i < nr_threads; i++)
        spm_mt->spm_threads[i].cpu = threads_cpus[i];

    ///> Load the appropriate sub-matrix to each thread.
    Spms = SPM::LoadMMF_mt(mmf_fname, nr_threads);

    ///> Start timer for preprocessing.
    xtimer_t timer;
    timer_init(&timer);
    timer_start(&timer);

    ///> Initalization of pthreads.
    threads = (pthread_t *) malloc((nr_threads - 1) * sizeof(pthread_t));
    if (!threads){
        perror("malloc");
        exit(1);
    }

    ///> Parameters passed to each thread.
    data = new Parameters[nr_threads];
    for (unsigned int i = 0; i < nr_threads; i++) {
        data[i].Spm = &Spms[i];
        data[i].wsize = wsize;
        data[i].thread_no = i;
        data[i].xform_buf = xform_buf;
        data[i].sampling_prob = sampling_prob;
        data[i].samples_max = samples_max;
        data[i].split_blocks = split_blocks;
        data[i].deltas = deltas;
        data[i].buffer.str("");
    }

    ///> Parallel Preprocessing for each thread.
    for (unsigned int i = 1; i < nr_threads; i++)
        pthread_create(&threads[i-1],NULL,thread_function,(void *) &data[i]);
        
    thread_function((void *) &data[0]);

    ///> Wait for other threads to finish.
    for (unsigned int i = 1; i < nr_threads; ++i)
        pthread_join(threads[i-1],NULL);

    CsxJitInitGlobal();

    ///> Init Jits.
    Jits = (CsxJit **) malloc(nr_threads * sizeof(CsxJit *));
    if (!Jits){
        perror("malloc");
        exit(1);
    }

    ///> Serial Preprocessing.
    for (unsigned int i = 0; i < nr_threads; ++i){
        /**
         *  Init CsxManager which holds the final CSX form of encoded
         *  sub-matrix.
         */
        CsxManager *CsxMg = new CsxManager(&Spms[i]);
        ///> Get the elements in CSX form.
        spm_mt->spm_threads[i].spm = CsxMg->MakeCsx();
        ///> Init Jits which makes the generation code part.
        Jits[i] = new CsxJit(CsxMg, i);
        ///> Make the generated code
        Jits[i]->genCode(data[i].buffer);
        //verifyModule(*(Jits[i]->M), AbortProcessAction, 0);
        delete CsxMg;
        ///> Print the results of each thread.
        std::cout << data[i].buffer.str();
    }

    ///> Optimize generated code and apply it to each thread seperately.
    CsxJitOptmize();
    for (unsigned int i=0; i < nr_threads; i++) {
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->getSpmvFn();
        delete Jits[i];
    }

    ///> Free and Delete unnecessary data.
    free(Jits);
    free(threads);
    delete[] Spms;
    delete[] data;

    ///> Stop the preprocessing time counter and print it.
    timer_pause(&timer);
    std::cout << "Preprocessing time: "
              << timer_secs(&timer) << " sec" 
              << std::endl;
    return spm_mt;
}

static void PutSpmMt(spm_mt_t *spm_mt)
{
    free(spm_mt->spm_threads);
    free(spm_mt);
}

///> Check if encoded matrix return the correct result.
static void CheckLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    void *crs;
    uint64_t nrows, ncols, nnz;

    crs = spm_crs32_double_init_mmf(mmf_name, &nrows, &ncols, &nnz);
    std::cout << "Checking ... " << std::flush;
    spmv_double_check_mt_loop(crs, spm_mt, spm_crs32_double_multiply, 1, nrows,
                              ncols, NULL);
    spm_crs32_double_destroy(crs);
    std::cout << "Check Passed" << std::endl << std::flush;
}

///> Count the size of the compressed matrix in CSX form.
static unsigned long CsxSize(spm_mt_t *spm_mt)
{
    unsigned long ret;

    ret = 0;
    for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        csx_double_t *csx = (csx_double_t *)t->spm;
        
        ret += csx->nnz * sizeof(double);
        ret += csx->ctl_size;
    }
    
    return ret;
}

///> Run spmv and print the results.
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
        spm_mt = GetSpmMt(argv[i]);
        CheckLoop(spm_mt, argv[i]);
        std::cerr.flush();
        BenchLoop(spm_mt, argv[i]);
        PutSpmMt(spm_mt);
    }
    return 0;
}


// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
