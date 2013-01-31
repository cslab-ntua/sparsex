/*
 * spmv.cc -- Front-end utilities for invoking CSX.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "spmv.h"


///> Max deltas that an encoding type may have. This is restricted by the
///  number of bits used for encoding the different patterns in CSX.
///
///  @see ctl_ll.h
#define DELTAS_MAX  CTL_PATTERNS_MAX

// malloc wrapper
#define xmalloc(x)                         \
({                                         \
    void *ret_;                            \
    ret_ = malloc(x);                      \
    if (ret_ == NULL){                     \
        std::cerr << __FUNCTION__          \
                  << " " << __FILE__       \
                  << ":" << __LINE__       \
                  << ": malloc failed\n";  \
        exit(1);                           \
    }                                      \
    ret_;                                  \
})

using namespace csx;

static int SplitString(char *str, char **str_buf, const char *start_sep,
                       const char *end_sep)
{
    char *token = strtok(str, start_sep);
    int next = 0;
    int str_length = strcspn(token, end_sep);

    str_buf[next] = (char *) xmalloc((str_length + 1) * sizeof(char));
    strncpy(str_buf[next], token, str_length);
    str_buf[next][str_length] = 0;
    ++next;
    while ((token = strtok(NULL, start_sep)) != NULL) {
        str_length = strcspn(token, end_sep);
        str_buf[next] = (char *) xmalloc((str_length - 1) * sizeof(char));
        strncpy(str_buf[next], token, str_length);
        str_buf[next][str_length] = 0;
        ++next;
    }

    return next;
}

static void GetOptionXform(int **xform_buf)
{
    char *xform_orig = getenv("XFORM_CONF");

    *xform_buf = (int *) xmalloc(XFORM_MAX * sizeof(int));
    if (xform_orig && strlen(xform_orig)) {
        int next;
        int t;
        char *token, *xform;

        // Copy environment variable to avoid being altered from strtok()
        xform = (char *) xmalloc(strlen(xform_orig)+1);
        strncpy(xform, xform_orig, strlen(xform_orig)+1);

        next = 0;
        t = atoi(strtok(xform, ","));
        (*xform_buf)[next] = t;
        ++next;
        while ((token = strtok(NULL, ",")) != NULL) {
            t = atoi(token);
            (*xform_buf)[next] = t;
            ++next;
        }

        (*xform_buf)[next] = -1;
        free(xform);
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

static void GetOptionEncodeDeltas(int ***deltas)
{
    char *encode_deltas_env = getenv("ENCODE_DELTAS");
    if (encode_deltas_env && strlen(encode_deltas_env)) {
        // Copy environment variable to avoid being altered from strtok()
        char *encode_deltas_str = (char *) xmalloc(strlen(encode_deltas_env)+1);
        strncpy(encode_deltas_str, encode_deltas_env,
                strlen(encode_deltas_env)+1);

        // Init matrix deltas.
        *deltas = (int **) xmalloc(XFORM_MAX * sizeof(int *));

        for (int i = 0; i < XFORM_MAX; i++) {
            (*deltas)[i] = (int *) xmalloc(DELTAS_MAX * sizeof(int));
        }

        // Fill deltas with the appropriate data.
        char **temp = (char **) xmalloc(XFORM_MAX * sizeof(char *));
        int temp_size = SplitString(encode_deltas_str, temp, "{", "}");

        for (int i = 0; i < temp_size; i++) {
            int j = 0;
            char *token = strtok(temp[i], ",");

            (*deltas)[i][j] = atoi(token);
            ++j;
            while ((token = strtok(NULL, ",")) != NULL) {
                (*deltas)[i][j] = atoi(token);
                ++j;
            }

            (*deltas)[i][j] = -1;
            free(temp[i]);
        }

        free(temp);

        // Print deltas
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
        free(encode_deltas_str);
    }
}

static uint64_t GetOptionWindowSize()
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

static uint64_t GetOptionSamples()
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

static int GetOptionOuterLoops()
{
    const char *loops_env = getenv("OUTER_LOOPS");
    int ret = 1;
    
    if (loops_env) {
        ret = atoi(loops_env);
        if (ret < 0)
            ret = 0;
    }
    
    return ret;
}

static double GetOptionPortion()
{
    const char *sampling_portion_str = getenv("SAMPLING_PORTION");
    double sampling_portion;

    if (!sampling_portion_str) {
        sampling_portion = 0.0;
        std::cout << "Sampling portion: Not set" << std::endl;
    }
    else {
        sampling_portion = atof(sampling_portion_str);
        std::cout << "Sampling portion: " << sampling_portion << std::endl;
    }

    return sampling_portion;
}

void MakeMap(spm_mt_t *spm_mt, SPMSym *spm_sym)
{
    spm_mt_thread *spm_thread;
    SPM *spm;
    unsigned int *count;
    bool **initial_map;
    map_t *map;
    uint32_t start, end, col, total_count, temp_count, limit;
    uint32_t ncpus = spm_mt->nr_threads;
    uint32_t n = spm_sym->GetLowerMatrix()->GetNrCols();
#ifdef SPM_NUMA
    int node;
#endif
    
    ///> Init initial_map.
    count = (unsigned int *) xmalloc(n * sizeof(unsigned int));
    initial_map = (bool **) xmalloc(ncpus * sizeof(bool *));
    for (unsigned int i = 0; i < ncpus; i++)
        initial_map[i] = (bool *) xmalloc(n * sizeof(bool));
    
    for (unsigned int i = 0; i < n; i++)
        count[i] = 0;
        
    for (unsigned int i = 0; i < ncpus; i++)
        for (unsigned int j = 0; j < n; j++)
            initial_map[i][j] = 0;

    ///> Fill initial map.
    for (unsigned int i = 0; i < ncpus; i++) {
        spm = spm_sym[i].GetLowerMatrix();

        start = spm->GetRowStart();
        end = spm->GetRowStart() + spm->GetNrRows();
        
        for (uint32_t j = 0; j < spm->GetNrRows(); j++) {
            for (SpmRowElem *elem = spm->RowBegin(j); elem != spm->RowEnd(j);
                 elem++) {
                col = elem->x;
                assert(col < end);
                if (col < start + 1 && !initial_map[i][col]) {
                    initial_map[i][col] = 1;
                    count[col]++;
                }
            }
        }
    }
    total_count = 0;
    for (unsigned int i = 0; i < n; i++)
        total_count += count[i];
    
    ///> Print initial map.
    /*
    for (unsigned int i = 0; i < ncpus; i++) {
        for (unsigned int j = 0; j < n; j++)
            std::cout << initial_map[i][j] << " ";
        std::cout << std::endl;
    }
    for (unsigned int i = 0; i < n; i++)
        std::cout << count[i] << " ";
    std::cout << std::endl;
    */
    
    ///> Make map.
    end = 0;
    for (unsigned int i = 0; i < ncpus - 1; i++) {
        spm_thread = spm_mt->spm_threads + i;
#ifdef SPM_NUMA
        node = spm_thread->node;
        map = (map_t *) numa_alloc_onnode(sizeof(map_t), node);
#else
        map = (map_t *) xmalloc(sizeof(map_t));
#endif

        start = end;
        limit = total_count / (ncpus - i);
        temp_count = 0;
        while (temp_count < limit)
            temp_count += count[end++];
        total_count -= temp_count;
        map->length = temp_count;
#ifdef SPM_NUMA
        map->cpus = (unsigned int *) 
                        numa_alloc_onnode(temp_count * sizeof(unsigned int),
                                          node);
        map->elems_pos = (unsigned int *)
                             numa_alloc_onnode(temp_count *
                                               sizeof(unsigned int), node);
#else
        map->cpus = (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
        map->elems_pos = 
            (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
#endif        
        temp_count = 0;
        for (unsigned int j = start; j < end; j++) {
            for (unsigned int k = 0; k < ncpus; k++) {
                if (initial_map[k][j] == 1) {
                     map->cpus[temp_count] = k;
                     map->elems_pos[temp_count++] = j - 1;
                }
            }
        }
        assert(temp_count == map->length);
        
        spm_thread->map = map;
    }
    start = end;
    end = n;
    temp_count = total_count;
    
    spm_thread = spm_mt->spm_threads + ncpus - 1;

#ifdef SPM_NUMA
    node = spm_thread->node;

    map = (map_t *) numa_alloc_onnode(sizeof(map_t), node);
    map->cpus = (unsigned int *)
                    numa_alloc_onnode(temp_count * sizeof(unsigned int), node);
    map->elems_pos = (unsigned int *)
                         numa_alloc_onnode(temp_count * sizeof(unsigned int),
                                           node);
#else
    map = (map_t *) xmalloc(sizeof(map_t));
    map->cpus = (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
    map->elems_pos = (unsigned int *)
                         xmalloc(temp_count * sizeof(unsigned int));
#endif

    map->length = temp_count;
    temp_count = 0;
    for (unsigned int j = start; j < end; j++) {
        for (unsigned int k = 0; k < ncpus; k++) {
            if (initial_map[k][j] == 1) {
                map->cpus[temp_count] = k;
                map->elems_pos[temp_count++] = j - 1;
            }
        }
    }
    assert(temp_count == map->length); 
    spm_thread->map = map;
#ifdef SPM_NUMA
    int alloc_err = 0;
    for (unsigned int i = 0; i < ncpus; i++) {
        node = spm_mt->spm_threads[i].node;
        map = spm_mt->spm_threads[i].map;
        
        alloc_err += check_region(map, sizeof(map_t), node);
        alloc_err += check_region(map->cpus, 
                                    map->length * sizeof(unsigned int), node);
        alloc_err += check_region(map->elems_pos,
                                    map->length * sizeof(unsigned int), node);
    }
    
    std::cout << "allocation check for map... "
              << ((alloc_err) ?
                 "FAILED (see above for more info)" : "DONE")
              << std::endl;

#endif
    ///> Print final map.
    /*
    for (unsigned int i = 0; i < ncpus; i++) {
        spm_thread = spm_mt->spm_threads + i;
        map = spm_thread->map;
        std::cout << "Thread " << i << std::endl;
        for (unsigned int j = 0; j < map->length; j++) {
            std::cout << "(" << map->cpus[j] << ", " << map->elems_pos[j]
                      << ")" << std::endl;
        }
        std::cout << std::endl;
    }
    */
    
    ///> Release initial map.
    for (unsigned int i = 0; i < ncpus; i++)
        free(initial_map[i]);
    free(initial_map);
    free(count);
}

uint64_t MapSize(void *spm)
{
    unsigned int i;
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread;
    uint64_t size = 0;
    
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        size += spm_thread->map->length * sizeof(uint32_t);
        size += spm_thread->map->length * sizeof(uint32_t);
        size += spm_thread->map->length * sizeof(double);
    }
    
    return size;
}

void *PreprocessThread(void *thread_info)
{
    thread_info_t *data = (thread_info_t *) thread_info;
    
    // Set CPU affinity.
    setaffinity_oncpu(data->cpu);
    
    // If symmetric option not set ...
    if (!data->symmetric) {
        DRLE_Manager *DrleMg;
        data->buffer << "==> Thread: #" << data->thread_no << std::endl;

        // Initialize the DRLE manager.
        DrleMg = new DRLE_Manager(data->spm, 4, 255, 0.05, data->wsize,
                                  DRLE_Manager::SPLIT_BY_NNZ,
                                  data->sampling_portion, data->samples_max,
                                  data->split_blocks, false);
                                  
        // Adjust the ignore settings properly.
        DrleMg->IgnoreAll();
        for (int i = 0; data->xform_buf[i] != -1; ++i)
            DrleMg->RemoveIgnore(static_cast<SpmIterOrder>(data->xform_buf[i]));

        // If the user supplies the deltas choices, encode the matrix
        // with the order given in XFORM_CONF, otherwise find
        // statistical data for the types in XFORM_CONF, choose the
        // best choise, encode it and proceed likewise until there is
        // no satisfying encoding.
        if (data->deltas)
            DrleMg->EncodeSerial(data->xform_buf, data->deltas);
        else
            DrleMg->EncodeAll(data->buffer);

        // DrleMg->MakeEncodeTree();

        csx_double_t *csx = data->csxmg->MakeCsx(false);
        data->spm_encoded->spm = csx;
        data->spm_encoded->nr_rows = csx->nrows;
        data->spm_encoded->row_start = csx->row_start;

#ifdef SPM_NUMA
        int alloc_err = 0;
        alloc_err = check_region(csx->ctl, csx->ctl_size * sizeof(uint8_t),
                                 data->spm_encoded->node);
        data->buffer << "allocation check for ctl field... "
                     << ((alloc_err) ?
                        "FAILED (see above for more info)" : "DONE")
                     << std::endl;

        alloc_err = check_region(csx->values, csx->nnz * sizeof(*csx->values),
                                 data->spm_encoded->node);
        data->buffer << "allocation check for values field... "
                     << ((alloc_err) ?
                        "FAILED (see above for more info)" : "DONE")
                     << std::endl;
#endif
        delete DrleMg;
    } else {  // If symmetric option is set split the sub-matrix into two.
        DRLE_Manager *DrleMg1;
        DRLE_Manager *DrleMg2;

        data->buffer << "==> Thread: #" << data->thread_no << std::endl;
        
        data->spm_sym->DivideMatrix();
        
        // Initialize the DRLE manager
        DrleMg1 = new DRLE_Manager(data->spm_sym->GetFirstMatrix(), 4, 255-1, 
                                   0.05, data->wsize,
                                   DRLE_Manager::SPLIT_BY_NNZ,
                                   data->sampling_portion, data->samples_max,
                                   false);
        DrleMg2 = new DRLE_Manager(data->spm_sym->GetSecondMatrix(), 4, 255-1, 
                                   0.05, data->wsize,
                                   DRLE_Manager::SPLIT_BY_NNZ,
                                   data->sampling_portion, data->samples_max,
                                   false);
                        
        // Adjust the ignore settings properly
        DrleMg1->IgnoreAll();
        DrleMg2->IgnoreAll();
        for (int i = 0; data->xform_buf[i] != -1; ++i) {
            DrleMg1->
                RemoveIgnore(static_cast<SpmIterOrder>(data->xform_buf[i]));
            DrleMg2->
                RemoveIgnore(static_cast<SpmIterOrder>(data->xform_buf[i]));
        }
        
        // If the user supplies the deltas choices, encode the matrix with the
        // order given in XFORM_CONF, otherwise find statistical data for the 
        // types in XFORM_CONF, choose the best choise, encode it and proceed 
        // likewise until there is no satisfying encoding.
        if (data->deltas) {
            DrleMg1->EncodeSerial(data->xform_buf, data->deltas);
            DrleMg2->EncodeSerial(data->xform_buf, data->deltas);
        } else {
            DrleMg1->EncodeAll(data->buffer);
            DrleMg2->EncodeAll(data->buffer);
        }
        data->spm_sym->MergeMatrix();

        csx_double_sym_t *csx = data->csxmg->MakeCsxSym();
        data->spm_encoded->spm = csx;
        data->spm_encoded->row_start = csx->lower_matrix->row_start;
        data->spm_encoded->nr_rows = csx->lower_matrix->nrows;

#ifdef SPM_NUMA
        int alloc_err;
        
        alloc_err = check_region(csx->lower_matrix->ctl,
                                 csx->lower_matrix->ctl_size * sizeof(uint8_t),
                                 data->spm_encoded->node);
        data->buffer << "allocation check for ctl field... "
                     << ((alloc_err) ?
                        "FAILED (see above for more info)" : "DONE")
                     << std::endl;

        alloc_err = check_region(csx->lower_matrix->values,
                                 csx->lower_matrix->nnz * sizeof(double),
                                 data->spm_encoded->node);
        data->buffer << "allocation check for values field... "
                     << ((alloc_err) ?
                        "FAILED (see above for more info)" : "DONE")
                     << std::endl;
                     
	alloc_err = check_region(csx->dvalues,
	                         csx->lower_matrix->nrows * sizeof(double),
	                         data->spm_encoded->node);
        data->buffer << "allocation check for dvalues field... "
                     << ((alloc_err) ?
                        "FAILED (see above for more info)" : "DONE")
                     << std::endl;
#endif

        delete DrleMg1;
        delete DrleMg2;
    }
    return 0;
}

spm_mt_t *GetSpmMt(char *mmf_fname, CsxExecutionEngine &engine,
                   bool split_blocks, bool symmetric, SPM *spms)
{
    unsigned int nr_threads, *threads_cpus;
    spm_mt_t *spm_mt;
    bool own_spms = !spms;
    int *xform_buf = NULL;
    int **deltas = NULL;
    thread_info_t *data;
    pthread_t *threads;
    SPMSym *spms_sym = NULL;

    // Get MT_CONF
    mt_get_options(&nr_threads, &threads_cpus);
    setaffinity_oncpu(threads_cpus[0]);
    std::cout << "MT_CONF=";
    for (unsigned int i = 0; i < nr_threads; i++) {
        if (i != 0)
            std::cout << ",";
        std::cout << threads_cpus[i];
    }
    std::cout << std::endl;

    // Get XFORM_CONF
    GetOptionXform(&xform_buf);

    // Get ENCODE_DELTAS
    GetOptionEncodeDeltas(&deltas);

    // Get WINDOW_SIZE
    uint64_t wsize = GetOptionWindowSize();

    // Get SAMPLES
    uint64_t samples_max = GetOptionSamples();

    // Get SAMPLING_PORTION
    double sampling_portion = GetOptionPortion();

    // Set affinity for the serial part of the preproprecessing.
    setaffinity_oncpu(threads_cpus[0]);

    // Initialization of the multithreaded sparse matrix representation
    spm_mt = (spm_mt_t *) xmalloc(sizeof(spm_mt_t));

    spm_mt->nr_threads = nr_threads;
    spm_mt->symmetric = symmetric;
    spm_mt->spm_threads =
        (spm_mt_thread_t *) xmalloc(sizeof(spm_mt_thread_t) * nr_threads);

    for (unsigned int i = 0; i < nr_threads; i++) {
        spm_mt->spm_threads[i].cpu = threads_cpus[i];
        spm_mt->spm_threads[i].node = numa_node_of_cpu(threads_cpus[i]);
        spm_mt->spm_threads[i].id = i;
    }

    // Load the appropriate sub-matrix to each thread
    if (!symmetric) {
        if (!spms)
            spms = SPM::LoadMMF_mt(mmf_fname, nr_threads);
    } else {
        if (!spms_sym) {
            spms_sym = SPMSym::LoadMMF_mt(mmf_fname, nr_threads);
            // Switch Reduction Phase
            MakeMap(spm_mt, spms_sym);
        }
    }
    
    // Start timer for preprocessing
    xtimer_t timer;
    timer_init(&timer);
    timer_start(&timer);

    // Initialize and setup threads
    threads = (pthread_t *) xmalloc(nr_threads * sizeof(pthread_t));

    data = new thread_info_t[nr_threads];
    for (unsigned int i = 0; i < nr_threads; i++) {
        if (!symmetric) {
            data[i].spm = spms + i;
            data[i].spm_sym = NULL;
            data[i].csxmg = new CsxManager(data[i].spm);
        } else {
            data[i].spm = NULL;
            data[i].spm_sym = spms_sym + i;
            data[i].csxmg = new CsxManager(spms_sym + i);
        }
        data[i].spm_encoded = &spm_mt->spm_threads[i];
        data[i].wsize = wsize;
        data[i].thread_no = i;
        data[i].cpu = threads_cpus[i];
        data[i].xform_buf = xform_buf;
        data[i].sampling_portion = sampling_portion;
        data[i].samples_max = samples_max;
        data[i].split_blocks = split_blocks;
        data[i].symmetric = symmetric;
        data[i].deltas = deltas;
        data[i].buffer.str("");
#ifdef SPM_NUMA
        // Enable the full-column-index optimization for NUMA architectures
        data[i].csxmg->SetFullColumnIndices(true);
#endif
    }
    
    // Start parallel preprocessing
    for (unsigned int i = 0; i < nr_threads; i++)
        pthread_create(&threads[i], NULL, PreprocessThread,
                       (void *) &data[i]);

    for (unsigned int i = 0; i < nr_threads; ++i)
        pthread_join(threads[i], NULL);

    // CSX JIT compilation
    CsxJit **Jits = new CsxJit*[nr_threads];
    for (unsigned int i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit(data[i].csxmg, &engine, i, symmetric);
        Jits[i]->GenCode(data[i].buffer);
        std::cout << data[i].buffer.str();
    }

    for (unsigned int i = 0; i < nr_threads; i++)
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();
    
    // Cleanup.
    free(threads);
    free(threads_cpus);
    free(xform_buf);
    
    delete[] Jits;
    delete[] data;
    if (own_spms)
        delete[] spms;

    // Preprocessing finished; stop timer
    timer_pause(&timer);
    pre_time = timer_secs(&timer);
    
    return spm_mt;
}

void PutSpmMt(spm_mt_t *spm_mt)
{
    if (!spm_mt->symmetric) {
        for (unsigned i = 0; i < spm_mt->nr_threads; i++) {
            csx_double_t *csx = (csx_double_t *) spm_mt->spm_threads[i].spm;
            
            DestroyCsx(csx);
        }
    } else {
        for (unsigned i = 0; i < spm_mt->nr_threads; i++) {
            csx_double_sym_t *csx_sym =
                (csx_double_sym_t *) spm_mt->spm_threads[i].spm;
            
            DestroyCsxSym(csx_sym);
        }
    }
    free(spm_mt->spm_threads);
    free(spm_mt);
}

void CheckLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    void *crs;
    uint64_t nrows, ncols, nnz;

    crs = spm_crs32_double_init_mmf(mmf_name, &nrows, &ncols, &nnz, NULL);
    std::cout << "Checking ... " << std::flush;
    if (!spm_mt->symmetric)
        SPMV_CHECK_FN(crs, spm_mt, spm_crs32_double_multiply, 1, nrows, ncols,
                      NULL);
    else
        SPMV_CHECK_SYM_FN(crs, spm_mt, spm_crs32_double_multiply, 1, nrows,
                          ncols, NULL);

    spm_crs32_double_destroy(crs);
    std::cout << "Check Passed" << std::endl;
}

uint64_t CsxSize(void *spm_mt)
{
    return (uint64_t) CsxSize((spm_mt_t *) spm_mt);
}

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

uint64_t CsxSymSize(void *spm_mt)
{
    return (uint64_t) CsxSymSize((spm_mt_t *) spm_mt);
}

static unsigned long CsxSymSize(spm_mt_t *spm_mt)
{
    unsigned long ret;

    ret = 0;
    for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t *t = spm_mt->spm_threads + i;
        csx_double_sym_t *csx_sym = (csx_double_sym_t *)t->spm;
        csx_double_t *csx = csx_sym->lower_matrix;
        
        ret += (csx->nnz + csx->nrows) * sizeof(double);
        ret += csx->ctl_size;
    }

    return ret;
}    

void BenchLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    uint64_t nrows, ncols, nnz;
    double secs, flops;
    long loops_nr = 128;

    ReadMmfSizeLine(mmf_name, nrows, ncols, nnz);
    int nr_outer_loops = GetOptionOuterLoops();
    
    for (int i = 0; i < nr_outer_loops; ++i) {
	if (!spm_mt->symmetric) {
            secs = SPMV_BENCH_FN(spm_mt, loops_nr, nrows, ncols, NULL);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            printf("m:%s f:%s s:%lu pt:%lf t:%lf r:%lf\n", "csx",
                   basename(mmf_name), CsxSize(spm_mt), pre_time, secs,
                   flops);
        } else {
            secs = SPMV_BENCH_SYM_FN(spm_mt, loops_nr, nrows, ncols, NULL);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            // Switch Reduction Phase
            printf("m:%s f:%s ms:%lu s:%lu pt:%lf t:%lf r:%lf\n", "csx-sym",
                   basename(mmf_name), MapSize(spm_mt), CsxSymSize(spm_mt),
                   pre_time, secs, flops);
        }        
    }
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
