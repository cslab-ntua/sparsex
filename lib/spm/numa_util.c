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
#include <assert.h>

void *alloc_interleaved(size_t size, size_t *parts, size_t nr_parts,
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

	struct bitmask *nodemask = numa_allocate_nodemask();
#define PAGE_ALIGN(addr) (void *)((unsigned long) addr & ~(pagesize-1))

	/*
	 * Bind parts to specific nodes
	 * All parts must be page aligned
	 */
	void *curr_part = ret;
	size_t i;
	size_t rem = 0;
	size_t new_part_size = 0;
	for (i = 0; i < nr_parts; i++) {
		size_t	part_size = parts[i] + rem;
		rem = part_size % pagesize;
		if (part_size < pagesize) {
			// FIXME: not stable when more than two partitions
			// fall in the same page
			new_part_size += part_size;
			if (new_part_size < pagesize) {
				parts[i] = 0;
				continue;
			} else {
				part_size = new_part_size;
			}
		} else {
			if (i < nr_parts - 1) {
				if (rem < pagesize / 2) {
					part_size -= rem;
				} else {
					part_size += pagesize - rem;
					// rem is used to adjust the size of
					// the next part
					rem -= pagesize;
				}
			}
		}

		numa_bitmask_setbit(nodemask, nodes[i]);
		if (mbind(PAGE_ALIGN(curr_part), part_size, MPOL_BIND,
		    nodemask->maskp, nodemask->size, 0) < 0) {
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

int check_interleaved(void *addr, size_t size, size_t *parts,
                       size_t nr_parts, const int *nodes)
{
	int node;
	unsigned long i;
	unsigned long cur_size = 0;
	int cur_part = 0;
	int pagesize = numa_pagesize();
	int ret = 1;
	void * start = NULL;
	void * end = NULL;
	int error_node = 0;
	int correct_node = 0;
	

	for (i = 0; i < size; i += pagesize, cur_size += pagesize) {
		if (get_mempolicy(&node, 0, 0, addr + i,
                                  MPOL_F_ADDR | MPOL_F_NODE) < 0) {
			perror("get_mempolicy()");
			exit(1);
		}

		if (cur_size == parts[cur_part]) {
			cur_size = 0;
			cur_part++;
		}

		if (node != nodes[cur_part]) {
			ret = 0;
			if (!start || end != addr + i - 1 ||
                            error_node != node ||
			    correct_node != nodes[cur_part]) {
				if (start) {
                        		fprintf(stderr,
					        "pages %p - %p are on node %d"
						" and must be on node %d\n",
					        start, end, error_node,
					        correct_node);
				}
				start = addr + i;
				end = addr + i + pagesize - 1;
				error_node = node;
				correct_node = nodes[cur_part];
			} else {
				end = addr + i + pagesize - 1;
			}
		}
	}
	if (start)
		fprintf(stderr,
		        "pages %p - %p are on node %d and must be on node %d\n",
		        start, end, error_node, correct_node);
	
	assert(cur_part == nr_parts - 1);
	assert(cur_size >= parts[cur_part] &&
	       cur_size - parts[cur_part] < pagesize);
	return ret;
}

int check_onnode(void *addr, size_t size, const int node)
{
	return check_interleaved(addr, size, &size, 1, &node);
}
