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

#if 0   // SYM
void MakeMap(spm_mt_t *spm_mt, SparsePartitionSym<uint64_t, double> *spm_sym)
{
    spm_mt_thread *spm_thread;
    SparsePartition<uint64_t, double> *spm;
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
            for (Elem<uint64_t, double> *elem = spm->RowBegin(j); elem != spm->RowEnd(j);
                 elem++) {
                col = elem->col;
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
#endif  // SYM    

spm_mt_t *PrepareSpmMt(const RuntimeContext &rt_config,
                       const CsxContext &csx_config)
{
    spm_mt_t *spm_mt;
    spm_mt = (spm_mt_t *) xmalloc(sizeof(spm_mt_t));

    spm_mt->nr_threads = rt_config.GetNrThreads();
    spm_mt->symmetric = csx_config.IsSymmetric();
    spm_mt->spm_threads = (spm_mt_thread_t *) xmalloc
        (sizeof(spm_mt_thread_t) * rt_config.GetNrThreads());

    for (size_t i = 0; i < rt_config.GetNrThreads(); i++) {
        spm_mt->spm_threads[i].cpu = rt_config.GetAffinity(i);
        spm_mt->spm_threads[i].node =
            numa_node_of_cpu(rt_config.GetAffinity(i));
        spm_mt->spm_threads[i].id = i;
    }

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

#if 0   // SYM
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
#endif  // SYM

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
        }/* else {
            secs = SPMV_BENCH_SYM_FN(spm_mt, loops_nr, nrows, ncols, NULL);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            // Switch Reduction Phase
            printf("m:%s f:%s ms:%lu s:%lu pt:%lf t:%lf r:%lf\n", "csx-sym",
                   basename(mmf_name), MapSize(spm_mt), CsxSymSize(spm_mt),
                   pre_time, secs, flops);
                   }
         */      
    }
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
