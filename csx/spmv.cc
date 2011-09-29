/*
 * spmv.cc -- Utility routines for invoking CSX.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2011, Vasileios Karakasis
 * Copyright (C) 2010-2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <numa.h>

#include "spm.h"
#include "mmf.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"
#include "llvm_jit_help.h"
#include "spmv.h"

#include <numa.h>
#include <numaif.h>

extern "C" {
#include "mt_lib.h"
#include "spmv_method.h"
#include "spm_crs.h"
#include "spm_mt.h"
#include "spmv_loops_mt.h"
#include "spmv_loops_sym_mt.h"
#ifdef SPM_NUMA
#   include "spmv_loops_mt_numa.h"
#else
#   include "spmv_loops_mt.h"
#endif
#include "timer.h"
#include "ctl_ll.h"
#include <libgen.h>
}

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

int SplitString(char *str, char **str_buf, const char *start_sep,
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

void GetOptionXform(int **xform_buf)
{
    char *xform_orig = getenv("XFORM_CONF");
    *xform_buf = (int *) xmalloc(XFORM_MAX * sizeof(int));

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

void GetOptionEncodeDeltas(int ***deltas)
{
    char *encode_deltas_str = getenv("ENCODE_DELTAS");

    if (encode_deltas_str) {
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
    }
}

uint64_t GetOptionWindowSize()
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

uint64_t GetOptionSamples()
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

double GetOptionProbability()
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

bool GetOptionSplitBlocks()
{
    const char *split_blocks_str = getenv("SPLIT_BLOCKS");
    bool split_blocks;

    if (!split_blocks_str)
        split_blocks = false;
    else
        split_blocks = true;

    return split_blocks;
}

bool GetOptionSymmetric()
{
    const char *symmetric_str = getenv("SYMMETRIC");
    bool symmetric;

    if (!symmetric_str)
        symmetric = false;
    else
        symmetric = true;

    return symmetric;
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
        map = (map_t *) xmalloc(sizeof(map_t));
        
        start = end;
        limit = total_count / (ncpus - i);
        temp_count = 0;
        while (temp_count < limit)
            temp_count += count[end++];
        total_count -= temp_count;
        /*std::cout << "Thread " << i << std::endl;
        std::cout << start << " -> " << end << std::endl;*/
        map->length = temp_count;
        map->cpus = (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
        map->elems_pos = 
            (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
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
    map = (map_t *) xmalloc(sizeof(map_t));
    
    /*std::cout << "Thread " << ncpus - 1 << std::endl;
    std::cout << start << " -> " << end << std::endl;*/
    
    map->length = temp_count;
    map->cpus = (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
    map->elems_pos = (unsigned int *) xmalloc(temp_count * sizeof(unsigned int));
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
    
    // set cpu affinity
    setaffinity_oncpu(data->cpu);
    if (!data->symmetric) {
        DRLE_Manager *DrleMg;
        data->buffer << "==> Thread: #" << data->thread_no << std::endl;

        // Initialize the DRLE manager
        DrleMg = new DRLE_Manager(data->spm, 4, 255-1, 0.1, data->wsize,
                                  DRLE_Manager::SPLIT_BY_NNZ,
                                  data->sampling_prob, data->samples_max);
                                  
        // Adjust the ignore settings properly
        DrleMg->IgnoreAll();
        for (int i = 0; data->xform_buf[i] != -1; ++i)
            DrleMg->RemoveIgnore(static_cast<SpmIterOrder>(data->xform_buf[i]));

        // If the user supplies the deltas choices, encode the matrix with the
        // order given in XFORM_CONF, otherwise find statistical data for the 
        // types in XFORM_CONF, choose the best choise, encode it and proceed 
        // likewise until there is no satisfying encoding.
        if (data->deltas)
            DrleMg->EncodeSerial(data->xform_buf, data->deltas,
                                 data->split_blocks);
        else
            DrleMg->EncodeAll(data->buffer, data->split_blocks);

        // DrleMg->MakeEncodeTree(data->split_blocks);

        csx_double_t *csx = data->csxmg->MakeCsx(false);
        data->spm_encoded->spm = csx;
        data->spm_encoded->nr_rows = csx->nrows;
        data->spm_encoded->row_start = csx->row_start;
#ifdef SPM_NUMA
        int node;
        if (get_mempolicy(&node, 0, 0, csx->values,
                          MPOL_F_ADDR | MPOL_F_NODE) < 0) {
            perror("get_mempolicy");
            exit(1);
        }

        data->buffer << "csx part " << data->thread_no
                     << " is on node " << node
                     << " and must be on node "
                     << data->spm_encoded->node << std::endl;

#endif
        delete DrleMg;
    } else {
        DRLE_Manager *DrleMg1;
        DRLE_Manager *DrleMg2;

        data->buffer << "==> Thread: #" << data->thread_no << std::endl;

        data->spm_sym->DivideMatrix();
        
        /*
        if (data->thread_no == 0) {
            std::cout << "Rows: " << data->spm_sym->GetFirstMatrix()->GetNrRows() << std::endl;
            std::cout << "Cols: " << data->spm_sym->GetFirstMatrix()->GetNrCols() << std::endl;
            std::cout << "Nnz: " << data->spm_sym->GetFirstMatrix()->GetNrNonzeros() << std::endl;
            std::cout << "Elems: " << data->spm_sym->GetFirstMatrix()->GetElemsSize() << std::endl;
            data->spm_sym->GetFirstMatrix()->PrintRows();
            data->spm_sym->GetFirstMatrix()->PrintElems();
            std::cout << "Rows: " << data->spm_sym->GetSecondMatrix()->GetNrRows() << std::endl;
            std::cout << "Cols: " << data->spm_sym->GetSecondMatrix()->GetNrCols() << std::endl;
            std::cout << "Nnz: " << data->spm_sym->GetSecondMatrix()->GetNrNonzeros() << std::endl;
            std::cout << "Elems: " << data->spm_sym->GetSecondMatrix()->GetElemsSize() << std::endl;
            data->spm_sym->GetSecondMatrix()->PrintRows();
            data->spm_sym->GetSecondMatrix()->PrintElems();
        }
        */
        
        // Initialize the DRLE manager
        DrleMg1 = new DRLE_Manager(data->spm_sym->GetFirstMatrix(), 4, 255-1, 
                                   0.1, data->wsize, DRLE_Manager::SPLIT_BY_NNZ,
                                   data->sampling_prob, data->samples_max);
        DrleMg2 = new DRLE_Manager(data->spm_sym->GetSecondMatrix(), 4, 255-1, 
                                   0.1, data->wsize, DRLE_Manager::SPLIT_BY_NNZ,
                                   data->sampling_prob, data->samples_max);
                                  
        // Adjust the ignore settings properly
        DrleMg1->IgnoreAll();
        DrleMg2->IgnoreAll();
        for (int i = 0; data->xform_buf[i] != -1; ++i) {
            DrleMg1->RemoveIgnore(static_cast<SpmIterOrder>(data->xform_buf[i]));
            DrleMg2->RemoveIgnore(static_cast<SpmIterOrder>(data->xform_buf[i]));
        }
        
        // If the user supplies the deltas choices, encode the matrix with the
        // order given in XFORM_CONF, otherwise find statistical data for the 
        // types in XFORM_CONF, choose the best choise, encode it and proceed 
        // likewise until there is no satisfying encoding.
        if (data->deltas) {
            DrleMg1->EncodeSerial(data->xform_buf, data->deltas,
                                  data->split_blocks);
            DrleMg2->EncodeSerial(data->xform_buf, data->deltas,
                                  data->split_blocks);
        } else {
            DrleMg1->EncodeAll(data->buffer, data->split_blocks);
            DrleMg2->EncodeAll(data->buffer, data->split_blocks);
        }
        
        /*
        if (data->thread_no == 0) {
            std::cout << "Rows: " << data->spm_sym->GetFirstMatrix()->GetNrRows() << std::endl;
            std::cout << "Cols: " << data->spm_sym->GetFirstMatrix()->GetNrCols() << std::endl;
            std::cout << "Nnz: " << data->spm_sym->GetFirstMatrix()->GetNrNonzeros() << std::endl;
            std::cout << "Elems: " << data->spm_sym->GetFirstMatrix()->GetElemsSize() << std::endl;
            data->spm_sym->GetFirstMatrix()->PrintRows();
            //data->spm_sym->GetFirstMatrix()->PrintElems();
            std::cout << "Rows: " << data->spm_sym->GetSecondMatrix()->GetNrRows() << std::endl;
            std::cout << "Cols: " << data->spm_sym->GetSecondMatrix()->GetNrCols() << std::endl;
            std::cout << "Nnz: " << data->spm_sym->GetSecondMatrix()->GetNrNonzeros() << std::endl;
            std::cout << "Elems: " << data->spm_sym->GetSecondMatrix()->GetElemsSize() << std::endl;
            data->spm_sym->GetSecondMatrix()->PrintRows();
            //data->spm_sym->GetSecondMatrix()->PrintElems();
        }
        */
        
        data->spm_sym->MergeMatrix();
        
        /*
        std::cout << "Rows: " << data->spm_sym->GetLowerMatrix()->GetNrRows() << std::endl;
        std::cout << "Cols: " << data->spm_sym->GetLowerMatrix()->GetNrCols() << std::endl;
        std::cout << "Nnz: " << data->spm_sym->GetLowerMatrix()->GetNrNonzeros() << std::endl;
        std::cout << "Elems: " << data->spm_sym->GetLowerMatrix()->GetElemsSize() << std::endl;
        data->spm_sym->GetLowerMatrix()->PrintRows();
        //data->spm_sym->GetLowerMatrix()->PrintElems();
        */
        
        csx_double_sym_t *csx = data->csxmg->MakeCsxSym();
        data->spm_encoded->spm = csx;
        data->spm_encoded->row_start = csx->lower_matrix->row_start;
        data->spm_encoded->nr_rows = csx->lower_matrix->nrows;
        
        delete DrleMg1;
        delete DrleMg2;
    }
    
    return 0;
}

spm_mt_t *GetSpmMt(char *mmf_fname, csx::SPM *spms)
{
    unsigned int nr_threads, *threads_cpus;
    spm_mt_t *spm_mt;
    int *xform_buf = NULL;
    int **deltas = NULL;
    thread_info_t *data;
    pthread_t *threads;
    CsxJit **Jits;
    SPMSym *spms_sym = NULL;

    // Get MT_CONF
    mt_get_options(&nr_threads, &threads_cpus);
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

    // Get SAMPLING_PROB
    double sampling_prob = GetOptionProbability();

    // Get SPLIT_BLOCKS
    bool split_blocks = GetOptionSplitBlocks();
    
    // Get SYMMETRIC
    bool symmetric = GetOptionSymmetric();

    // Initalization of the multithreaded sparse matrix representation
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
            MakeMap(spm_mt, spms_sym);
        }
    }

    // Start timer for preprocessing
    xtimer_t timer;
    timer_init(&timer);
    timer_start(&timer);

    // Initalize and setup threads
    threads = (pthread_t *) xmalloc((nr_threads - 1) * sizeof(pthread_t));

    data = new thread_info_t[nr_threads];
    for (unsigned int i = 0; i < nr_threads; i++) {
        if (!symmetric) { 
            data[i].spm = spms + i;
            data[i].spm_sym = NULL;
            data[i].csxmg = new CsxManager(data[i].spm);
        }
        else {
            data[i].spm = NULL;
            data[i].spm_sym = spms_sym + i;
            data[i].csxmg = new CsxManager(spms_sym + i);
        }
        data[i].spm_encoded = &spm_mt->spm_threads[i];
        data[i].wsize = wsize;
        data[i].thread_no = i;
        data[i].cpu = threads_cpus[i];
        data[i].xform_buf = xform_buf;
        data[i].sampling_prob = sampling_prob;
        data[i].samples_max = samples_max;
        data[i].split_blocks = split_blocks;
        data[i].symmetric = symmetric;
        data[i].deltas = deltas;
        data[i].buffer.str("");
    }

    // Start parallel preprocessing
    for (unsigned int i = 1; i < nr_threads; i++)
        pthread_create(&threads[i-1], NULL, PreprocessThread,
                       (void *) &data[i]);

    PreprocessThread((void *) &data[0]);

    for (unsigned int i = 1; i < nr_threads; ++i)
        pthread_join(threads[i-1],NULL);

    /*
    for (unsigned int i = 0; i < nr_threads; ++i) {
        if (!symmetric) {
            csx_double_t *csx = (csx_double_t *) data[i].spm_encoded->spm;
            
            std::cout << "Rows: " << csx->nrows << std::endl;
            std::cout << "Cols: " << csx->ncols << std::endl;
            std::cout << "Non-Zero Elements: " << csx->nnz << std::endl;
            std::cout << "Row Start: " << csx->row_start << std::endl;
            std::cout << "Ctl Size: " << csx->ctl_size << std::endl;
            for (uint64_t i = 0; i < csx->ctl_size; i++)
                std::cout << (uint32_t) csx->ctl[i] << std::endl;
            
        } else {
            csx_double_sym_t *csx_sym = (csx_double_sym_t *) data[i].spm_encoded->spm;
            csx_double_t *csx = csx_sym->lower_matrix;
            
            std::cout << "Rows: " << csx->nrows << std::endl;
            std::cout << "Cols: " << csx->ncols << std::endl;
            std::cout << "Non-Zero Elements: " << csx->nnz << std::endl;
            std::cout << "Row Start: " << csx->row_start << std::endl;
            std::cout << "Ctl Size: " << csx->ctl_size << std::endl;
            std::cout << "Ctl: " << std::endl;
            for (uint64_t i = 0; i < csx->ctl_size; i++)
                std::cout << (uint32_t) csx->ctl[i] << std::endl;
            std::cout << "Values: " << std::endl;
            for (uint64_t i = 0; i < csx->nnz; i++)
                std::cout << csx->values[i] << std::endl;
            std::cout << "Diagonal Values: " << std::endl;
            for (uint64_t i = 0; i < csx->nrows; i++)
                std::cout << csx_sym->dvalues[i] << std::endl;
            
        }
    }
    */
        
    // CSX matrix construction and JIT compilation
    CsxJitInitGlobal();
    Jits = (CsxJit **) xmalloc(nr_threads * sizeof(CsxJit *));

    for (unsigned int i = 0; i < nr_threads; ++i) {
        Jits[i] = new CsxJit(data[i].csxmg, i, symmetric);
        Jits[i]->GenCode(data[i].buffer, symmetric);
        std::cout << data[i].buffer.str();
    }

    // Optimize generated code and assign it to every thread
    CsxJitOptmize();
    for (unsigned int i = 0; i < nr_threads; i++) {
        spm_mt->spm_threads[i].spmv_fn = Jits[i]->GetSpmvFn();
        delete Jits[i];
    }
    
    // Cleanup.
    free(Jits);
    free(threads);
    free(threads_cpus);
    free(xform_buf);
    delete[] spms;
    delete[] data;

    // Preprocessing finished; stop timer
    timer_pause(&timer);
    std::cout << "Preprocessing time: "
              << timer_secs(&timer) << " sec"
              << std::endl;
    return spm_mt;
}

void PutSpmMt(spm_mt_t *spm_mt)
{
    free(spm_mt->spm_threads);
    free(spm_mt);
}

void CheckLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    void *crs;
    uint64_t nrows, ncols, nnz;

    crs = spm_crs32_double_init_mmf(mmf_name, &nrows, &ncols, &nnz);
    std::cout << "Checking ... " << std::flush;

#ifdef SPM_NUMA
    spmv_double_check_mt_loop_numa(crs, spm_mt, spm_crs32_double_multiply, 2,
                                   nrows, ncols, NULL);
#else
    if (!spm_mt->symmetric)
        spmv_double_check_mt_loop(crs, spm_mt, spm_crs32_double_multiply, 2,
                                  nrows, ncols, NULL);
    else
        spmv_double_check_sym_mt_loop(crs, spm_mt, spm_crs32_double_multiply, 2,
                                      nrows, ncols, NULL);
#endif
    spm_crs32_double_destroy(crs);
    std::cout << "Check Passed" << std::endl << std::flush;
}

uint64_t CsxSize(void *spm_mt)
{
    return (uint64_t) CsxSize((spm_mt_t *) spm_mt);
}

unsigned long CsxSize(spm_mt_t *spm_mt)
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

unsigned long CsxSymSize(spm_mt_t *spm_mt)
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

    getMmfHeader(mmf_name, nrows, ncols, nnz);

#ifdef SPM_NUMA
    secs = spmv_double_bench_mt_loop_numa(spm_mt, loops_nr, nrows, ncols, NULL);
    flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
    printf("m:%s f:%s s:%lu t:%lf r:%lf\n",
           "csx", basename(mmf_name), CsxSize(spm_mt), secs, flops);
#else
    if (!spm_mt->symmetric) {
        secs = spmv_double_bench_mt_loop(spm_mt, loops_nr, nrows, ncols, NULL);
        flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
        printf("m:%s f:%s s:%lu t:%lf r:%lf\n",
               "csx", basename(mmf_name), CsxSize(spm_mt), secs, flops);
    } else {
        std::cout << "Map Size " << MapSize(spm_mt) << std::endl;
        secs = spmv_double_bench_sym_mt_loop(spm_mt, loops_nr, nrows, ncols,
                                             NULL);
        flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
        printf("m:%s f:%s s:%lu t:%lf r:%lf\n",
               "csx", basename(mmf_name), CsxSymSize(spm_mt), secs, flops);
    }
#endif
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
