/*
 * numa_util.c -- NUMA utilitiy functions.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
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

#define ALIGN_ADDR(addr, bound) (void *)((unsigned long) addr & ~(bound-1))

static int get_best_node(const size_t *node_scores, int nr_nodes)
{
    int j;
    int ret = 0;
    size_t max_score = node_scores[0];
    for (j = 1; j < nr_nodes; j++)
        if (node_scores[j] > max_score) {
            max_score = node_scores[j];
            ret = j;
        }

    return ret;
}

/**
 *  Fixes the user-specified interleaving to meet system's
 *  requirements. The updated partition sizes and allocation nodes are
 *  returned in `parts' and `nodes', respectively.
 *
 *  The adjustment is performed in two phases:
 *
 *  (a) Merging of very small partitions (less than the system's page
 *      size), in order to meet the system's allocation
 *      granularity. The partitions are merged from left to right,
 *      starting from a small partition and continuing until the size
 *      of the formed partition exceeds the system's page size. The
 *      `parts' array is updated accordingly; zeros are set in merged
 *      partitions. Each newly formed partition is assigned at the
 *      node with the largest allocation size according to the
 *      original partition list.
 *  (b) Fixing the partition boundaries, so that each partition starts
 *      at a multiple of the system's page size.
 */ 
static void fix_interleaving(size_t nr_parts, size_t *parts, int *nodes)
{
    size_t i;
	int nr_nodes = numa_num_configured_nodes();
    int pagesize = numa_pagesize();
    size_t *node_scores = calloc(nr_nodes, sizeof(*node_scores));
    if (!node_scores) {
        perror("malloc");
        exit(1);
    }

    // merge small partitions
    size_t base_part = 0;
    node_scores[nodes[0]] = parts[0];
    for (i = 1; i < nr_parts; i++) {
        size_t part_score = parts[i];
        if (parts[base_part] < pagesize) {
            // merge with base partition
            parts[base_part] += parts[i];
            parts[i] = 0;
        } else {
            // assign partition to the node with the highest score
            nodes[base_part] = get_best_node(node_scores, nr_nodes);
            
            // reset the scores
            memset(node_scores, 0, nr_nodes*sizeof(*node_scores));
            base_part = i;
        }

        node_scores[nodes[i]] += part_score;
    }

    // fix the last merger
    nodes[base_part] = get_best_node(node_scores, nr_nodes);

    // adjust partition size to multiples of system's page size
    size_t rem = 0;
    ssize_t size_adjust = 0;   // can be negative
    for (i = 0; i < nr_parts; i++) {
        if (!parts[i])
            continue;

        parts[i] += size_adjust;   // adjustment from previous partition
        rem = parts[i] % pagesize;

        // find next valid partition
        size_t next_part = i+1;
        while (next_part != nr_parts && !parts[next_part])
            next_part++;

        //
        // If at last partition, always execute the else path.
        // Otherwise, take the remainder of the page, if the largest part is
        // in the current partition, assuring, however, that the next
        // partition will have at least one page.
        // 
        if (next_part != nr_parts &&
            ((rem < pagesize / 2) ||
             (parts[next_part] + rem < 2*pagesize))) {
            parts[i] -= rem;
            size_adjust = rem;
        } else {
            parts[i] += pagesize - rem;
            size_adjust = rem - pagesize;
        }
    }

    free(node_scores);
}

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
void *alloc_interleaved(size_t size, size_t *parts, size_t nr_parts, int *nodes)
{
	size_t i;
	void *ret = NULL;

	/* sanity check */
	if (numa_available() < 0) {
		fprintf(stderr, "numa not implemented\n");
		exit(1);
	}

	ret = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE,
	           0, 0);
	if (ret == (void *) -1) {
		ret = NULL;
		goto exit;
	}

	struct bitmask *nodemask = numa_allocate_nodemask();

	// Fix the interleaving and bind parts to specific nodes
	fix_interleaving(nr_parts, parts, nodes);

	void *curr_part = ret;
	for (i = 0; i < nr_parts; i++) {
        if (!parts[i])
            continue;

		// Bind part to the proper node.
		numa_bitmask_setbit(nodemask, nodes[i]);
		if (parts[i] != 0 &&
            mbind(curr_part, parts[i], MPOL_BIND,
                  nodemask->maskp, nodemask->size, 0) < 0) {
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
		if (check_region(addr, parts[i], nodes[i]))
			err = 1;

		addr += parts[i];
	}

	return err;
}

void print_alloc_status(const char *data_descr, int err)
{
	printf("allocation check for %s... %s\n", data_descr,
	       (err) ? "FAILED (see above for more info)" : "DONE");
}
