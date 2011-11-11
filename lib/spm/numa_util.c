/*
 * numa_util.c -- NUMA utilitiy functions
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "numa_util.h"
#include <numa.h>
#include <numaif.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>

void *alloc_interleaved(size_t size,
                        size_t *parts, size_t nr_parts,
                        const int *nodes)
{
    void *ret = NULL;

    /* sanity check */
    if (numa_available() < 0) {
        fprintf(stderr, "numa not implemented\n");
        exit(1);
    }
    
	int pagesize = numa_pagesize();

    ret = mmap(NULL, size, PROT_READ | PROT_WRITE,
               MAP_ANONYMOUS | MAP_PRIVATE, 0, 0);
    if (ret == (void *) -1) {
        ret = NULL;
        goto exit;
    }

	struct bitmask *nodemask = numa_bitmask_alloc(numa_num_configured_cpus());
#define PAGE_ALIGN(addr) (void *)((unsigned long) addr & ~(pagesize-1))
	/*
	 * Bind parts to specific nodes
	 * All parts must be page aligned
	 */
	void *curr_part = ret;
    size_t i;
	size_t new_part_size = 0;
	for (i = 0; i < nr_parts; i++) {
		size_t	part_size = parts[i];
		size_t	rem = part_size % pagesize;
		if (part_size < pagesize) {
			new_part_size += part_size;
			if (new_part_size < pagesize) {
				parts[i] = 0;
				continue;
			} else {
				part_size = new_part_size;
			}
		} else {
			while (rem && rem < pagesize / 2 && i < nr_parts - 1) {
				/* Leave the page for the next partition */
				part_size--;
				rem = part_size % pagesize;
			}
		}

		numa_bitmask_setbit(nodemask, nodes[i]);
		if (mbind(PAGE_ALIGN(curr_part), part_size,
			      MPOL_BIND, nodemask->maskp, nodemask->size, 0) < 0) {
			perror("mbind");
			exit(1);
		}

		/* Clear the mask for the next round */
		numa_bitmask_clearbit(nodemask, nodes[i]);
		parts[i] = part_size;
		curr_part += parts[i];
		new_part_size = 0;
	}

#undef PAGE_ALIGN
	numa_bitmask_free(nodemask);
    
exit:
    return ret;
}

void free_interleaved(void *addr, size_t length)
{
    if (munmap(addr, length) < 0) {
        perror("munmap");
        exit(1);
    }
}
