/*
 * numa_util.c -- NUMA utilitiy functions
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * Copyright (C) 2011, Theodoros Gkountouvas
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

#define ALIGN_ADDR(addr, bound) (void *)((unsigned long) addr & ~(bound-1))

/**
 *  Allocate contiguous memory with a custom interleaving of pages on
 *  physical nodes.
 *
 *  @param size the total size of the allocation (in bytes). It will
 *              be rounded to the system's page size.
 *
 *  @param parts the size (in bytes) of each interleaving
 *               partition. The size of each partition will be
 *               automatically set to the nearest multiple of the page
 *               size. On return, parts are updated with the new
 *               sizes.
 *  @param nr_parts the number of interleaving parts.
 *  @param nodes the physical memory nodes to bind each partition. The
 *               size of the array must equal nr_parts.
 *  @return On success a pointer to the newly allocated area is
 *          returned, otherwise NULL is returned. 
 */
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
	if (ret == (void *) -1) {
		ret = NULL;
		goto exit;
	}

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
exit:
	return ret;
}

void *alloc_onnode(size_t size, int node)
{
    return alloc_interleaved(size, &size, 1, &node);
}

void free_interleaved(void *addr, size_t length)
{
	if (munmap(addr, length) < 0) {
		perror("munmap");
		exit(1);
	}
}

int check_region(void *addr, size_t size, int node)
{
	void *misplaced_start;
	int pagesize;
	size_t i;
	int misplaced_node;
	int err = 0;

	pagesize = numa_pagesize();
	misplaced_start = NULL;
	misplaced_node = -1;
	void *aligned_addr = ALIGN_ADDR(addr, pagesize);
	for (i = 0; i < size; i += pagesize) {
		int page_node;
		if (get_mempolicy(&page_node, 0, 0, aligned_addr + i,
		                  MPOL_F_ADDR | MPOL_F_NODE) < 0) {
			perror("get_mempolicy()");
			exit(1);
		}

		if (page_node != node) {
			// Start of a misplaced region
			if (!misplaced_start)
				misplaced_start = aligned_addr + i;
			misplaced_node = page_node;
			err = 1;
		} else {
			if (misplaced_start) {
				// End of a misplaced region
				assert(misplaced_node != -1);
				size_t misplaced_size = (aligned_addr + i - misplaced_start);
				fprintf(stderr, "Region [%p,%p) (%zd bytes) is misplaced "
				        "(lies on node %d but it should be on node %d)\n",
				        misplaced_start, aligned_addr + i, misplaced_size,
				        misplaced_node, node);
				misplaced_start = NULL;
				misplaced_node = -1;
			}
		}
	}

	if (misplaced_start) {
		// Last misplaced region
		assert(misplaced_node != -1);
		size_t misplaced_size = (aligned_addr + i - misplaced_start);
		fprintf(stderr, "Region [%p,%p) (%zd bytes) is misplaced "
		        "(lies on node %d but it should be on node %d)\n",
		        misplaced_start, aligned_addr + i, misplaced_size,
		        misplaced_node, node);
	}

	return err;
}

int check_interleaved(void *addr, const size_t *parts, size_t nr_parts,
                      const int *nodes)
{
	assert(addr && "addr is NULL");
	assert(parts && "parts is NULL");

	size_t i = 0;
	int err = 0;
	for (i = 0; i < nr_parts; i++) {
		if (check_region(addr, parts[i], nodes[i])) {
			err = 1;
		}
		addr += parts[i];
	}

	return err;
}

void print_alloc_status(const char *data_descr, int err)
{
	printf("allocation check for %s... %s\n", data_descr,
	       (err) ? "FAILED (see above for more info)" : "DONE");
}
