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

static void fix_interleaving(size_t nr_parts, size_t *parts, const int *nodes)
{
	int i, rem;
	int curr_page_size, max_page_size;
	int curr_part, chosen_node;
	int pagesize = numa_pagesize();
	int nr_nodes = numa_num_possible_nodes();
	int *page_to_nodes = (int *) malloc(nr_nodes * sizeof(int));

	for (curr_part = 0; curr_part < nr_parts - 1; curr_part++) {
		rem = parts[curr_part] % pagesize;

		// If a page is divided to more than one cpu ...
		if (rem != 0) {
                        // Initialize parameters.
                        for (i = 0; i < nr_nodes; i++)
                                page_to_nodes[i] = 0;
                        page_to_nodes[nodes[curr_part]] = rem;
                        curr_page_size = rem;
                        parts[curr_part] -= rem;

                        // Calculate page partitions to nodes.
			for (i = 1; curr_part + i < nr_parts &&
			     curr_page_size + parts[curr_part+i] < pagesize;
			     i++) {
				page_to_nodes[nodes[curr_part+i]] +=
				    parts[curr_part+i];
                                curr_page_size += parts[curr_part+i];
                                parts[curr_part+i] = 0;
                        }
                        if (curr_part + i < nr_parts) {
                                page_to_nodes[nodes[curr_part+i]] +=
				    pagesize - curr_page_size;
                                parts[curr_part+i] -=
				    pagesize - curr_page_size;
                                curr_page_size = pagesize;
                        }

                        // Assign current page to the proper node.
                        max_page_size = page_to_nodes[0];
                        chosen_node = 0;
                        for (i = 1; i < nr_nodes; i++) {
                                if (max_page_size < page_to_nodes[i]) {
                                        max_page_size = page_to_nodes[i];
                                        chosen_node = i;
                                }
                        }

                        // Assign current page to the proper CPU.
			for (i = 0; nodes[curr_part+i] != chosen_node; i++);
                        parts[curr_part+i] += curr_page_size;
		}
	}
	
}

void *alloc_interleaved(size_t size, size_t *parts, size_t nr_parts,
                        const int *nodes)
{
	int i;
	void *ret = NULL;

	/* sanity check */
	if (numa_available() < 0) {
		fprintf(stderr, "numa not implemented\n");
		exit(1);
	}

	ret = mmap(NULL, size, PROT_READ | PROT_WRITE,
                   MAP_ANONYMOUS | MAP_PRIVATE, 0, 0);

	if (ret == (void *) -1)
		return NULL;

	struct bitmask *nodemask = numa_allocate_nodemask();



	/*
	 * Bind parts to specific nodes
	 * All parts must be page aligned
	 */
	fix_interleaving(nr_parts, parts, nodes);

        void *curr_part = ret;
	for (i = 0; i < nr_parts; i++) {
		// Bind part to the proper node.
		numa_bitmask_setbit(nodemask, nodes[i]);
		
		if (parts[i] != 0 &&
		    mbind(curr_part, parts[i], MPOL_BIND, nodemask->maskp,
		          nodemask->size, 0) < 0) {
			perror("mbind");
			exit(1);
		}
		curr_part += parts[i];
		
		/* Clear the mask for the next round */
		numa_bitmask_clearbit(nodemask, nodes[i]);
	}
	numa_bitmask_free(nodemask);
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

		while (cur_size == parts[cur_part]) {
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
	
	assert(cur_size >= parts[cur_part] &&
	       cur_size - parts[cur_part] < pagesize);
	return ret;
}

int check_onnode(void *addr, size_t size, const int node)
{
	return check_interleaved(addr, size, &size, 1, &node);
}
