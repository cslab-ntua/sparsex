/*
 *	blocks.c -- Implementation of block operations and the block cache.
 *
 *	V. K. Karakasis, 2008, bkk@cslab.ece.ntua.gr
 */
#include "blocks.h"

#include <assert.h>
#include <stdlib.h>
#include "util.h"

/**
 *	List of block layouts and basic constants.
 */
struct storage_layout_s blk_layout_list[] = {
	{
		.sl_r = 0,
		.sl_c = 0,
		.sl_layout = BLOCK_STORAGE_RW,
		.sl_simd = 0
	},   /* The default layout */

	{
		.sl_r = 2,
		.sl_c = 2,
		.sl_layout = BLOCK_STORAGE_CW,
		.sl_simd = 1
	},

	{
		.sl_r = 2,
		.sl_c = 3,
		.sl_layout = BLOCK_STORAGE_CW,
		.sl_simd = 1
	},

	{
		.sl_r = 2,
		.sl_c = 4,
#if ELEM_TYPE == double
		.sl_layout = BLOCK_STORAGE_CW,
#else
		.sl_layout = BLOCK_STORAGE_RW,
#endif
		.sl_simd = 1
	},

	{
		.sl_r = 4,
		.sl_c = 2,
		.sl_layout = BLOCK_STORAGE_CW,
		.sl_simd = 1
	}
};

#define BLK_LAYOUT_LIST_SIZE                           \
    (sizeof(blk_layout_list) / sizeof(*blk_layout_list))

const int blk_layout_default = 0;
const int blk_layout_list_size = BLK_LAYOUT_LIST_SIZE;
const int blk_layout_list_begin = 0;
const int blk_layout_list_end = BLK_LAYOUT_LIST_SIZE - 1;

#if SPM_CRS_IDX_TYPE == int
/* `SPM_CRS_IDX_TYPE' is an `int', so use the fast compare function. */

static int
spm_index_cmp(const void *a, const void *b)
{
	SPM_CRS_IDX_TYPE *_a = (SPM_CRS_IDX_TYPE *) a;
	SPM_CRS_IDX_TYPE *_b = (SPM_CRS_IDX_TYPE *) b;

	/* This is a safe cast. We know `SPM_CRS_IDX_TYPE' is an `int'. */
	return (int) (*_a - *_b);
}

#else	/* SPM_CRS_IDX_TYPE != int */

static int
spm_index_cmp(const void *a, const void *b)
{
	SPM_CRS_IDX_TYPE *_a = (SPM_CRS_IDX_TYPE *) a;
	SPM_CRS_IDX_TYPE *_b = (SPM_CRS_IDX_TYPE *) b;

	if (*_a < *_b)
		return -1;
	else if (*_a == *_b)
		return 0;
	else
		return 1;
}

#endif	/* SPM_INDEX_BINDING */

static SPM_CRS_IDX_TYPE
map_cw(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c,
       SPM_CRS_IDX_TYPE i, SPM_CRS_IDX_TYPE j)
{
	return j*r + i;
}

static SPM_CRS_IDX_TYPE
map_rw(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c,
       SPM_CRS_IDX_TYPE i, SPM_CRS_IDX_TYPE j)
{
	return i*c + j;
}

block_t *
block_create(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c, int own)
{
	block_t	 *new = malloc(sizeof(*new));

	new->r   = r;
	new->c   = c;
	new->row = SPM_INDEX_INVALID;
	new->col = SPM_INDEX_INVALID;
	new->pad = r*c; /* the block is full of zeros */
	new->own = own;

	if (own)
		/* Allocate your own values. */
		new->values = calloc(r*c, sizeof(ELEM_TYPE));
	else
		new->values = NULL;

	return new;
}

void
block_init(block_t *block, SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c)
{
	block->r = r;
	block->c = c;
	block->row = SPM_INDEX_INVALID;
	block->col = SPM_INDEX_INVALID;
	block->pad = r*c;
	block->own = 0;
}

SPM_CRS_IDX_TYPE
block_get_nr_zeros(const block_t *block)
{
	SPM_CRS_IDX_TYPE nr_zeros;
	SPM_CRS_IDX_TYPE i;

	for (i = 0, nr_zeros = 0; i < block->r * block->c; i++)
		if (block->values[i] == 0.0)
			nr_zeros++;

	return (nr_zeros - block->pad);
}

void
block_delete(block_t *block)
{
	if (block && block->own)
		free(block->values);

	free(block);
	return;
}

#if 0
/*
 *	FIXME: Should use the `SPM_INDEX_FMT' and `SPM_VALUE_FMT' constants
 *	for the output format.
 */
void
block_dump(const block_t *block)
{
	SPM_CRS_IDX_TYPE i;

	fprintf(stderr, "[r = %" SPM_INDEX_FMT ", c = %" SPM_INDEX_FMT
	        ", row = %" SPM_INDEX_FMT ", col = %" SPM_INDEX_FMT
	        ", pad = %" SPM_INDEX_FMT ", ",
	        block->r, block->c, block->row, block->col, block->pad);

	fprintf(stderr, "values = { ");
	for (i = 0; i < block->r * block->c; i++)
		fprintf(stderr, "%.4lf ", block->values[i]);

	fprintf(stderr, "}]\n");
	return;
}
#endif

/*
 *	Block cache operations.
 */
block_cache_t *
block_cache_create(SPM_CRS_IDX_TYPE size,
                   SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c, int storage)
{
	block_cache_t	*new;

	new = malloc(sizeof(*new));
	new->blocks		= malloc(size*sizeof(*new->blocks));
	/* values are initialized to zero */
	new->blk_values = calloc(size*r*c, sizeof(ELEM_TYPE));
	new->valid		= bs_create(size);
	new->valid_pos	= calloc(size, sizeof(SPM_CRS_IDX_TYPE));
	new->r			= r;
	new->c			= c;
	new->size		= size;
	new->nr_values	= r*c*size;
	new->nr_blocks	= 0;

	/* Set up the blocks */
	SPM_CRS_IDX_TYPE i;
	for (i = 0; i < size; i++) {
		block_init(&new->blocks[i], r, c);
		/* Set the block values */
		new->blocks[i].values = &new->blk_values[i*r*c];
	}

	if (storage == BLOCK_STORAGE_CW)
		new->blk_map = map_cw;
	else
		new->blk_map = map_rw;

	return new;
}

void
block_cache_delete(block_cache_t *cache)
{
	if (!cache)
		return;

	bs_delete(cache->valid);
	free(cache->valid_pos);
	free(cache->blocks);
	free(cache->blk_values);
	free(cache);
	return;
}

void
block_cache_sort(block_cache_t *cache)
{
	qsort(cache->valid_pos, cache->nr_blocks, sizeof(SPM_CRS_IDX_TYPE),
		  spm_index_cmp);

	return;
}

void
block_cache_clear(block_cache_t *cache)
{
	SPM_CRS_IDX_TYPE i;

	for (i = 0; i < cache->nr_blocks; i++)
		bs_clear(cache->valid, cache->valid_pos[i]);

	cache->nr_blocks = 0;
	return;
}

int
block_cache_add_elem(block_cache_t *cache, const spm_elem_t *elem)
{
	SPM_CRS_IDX_TYPE bcol;
	SPM_CRS_IDX_TYPE brow;
	SPM_CRS_IDX_TYPE br, bc;
	int			is_new;

	bcol = elem->col / cache->c;
	bc   = elem->col % cache->c;
	brow = elem->row / cache->r;
	br   = elem->row % cache->r;

	is_new = !bs_isset(cache->valid, bcol);

	if (is_new) {
		block_clear(&cache->blocks[bcol]);
		bs_set(cache->valid, bcol);
		cache->valid_pos[cache->nr_blocks++] = bcol;
		cache->blocks[bcol].row = brow;
		cache->blocks[bcol].col = bcol;
	}

	cache->blk_values[bcol*cache->c*cache->r +
	                  cache->blk_map(cache->r, cache->c, br, bc)] = elem->value;

#if 0
	if (cache->blocks[bcol].pad <= 0) {
		DEBUG_PRINT_LONG(bcol);
		block_dump(&cache->blocks[bcol]);
	}
#endif

	assert(cache->blocks[bcol].pad);
	if (cache->blocks[bcol].pad != 0)
		cache->blocks[bcol].pad--;

	return is_new;
}

void
csr_add(SPM_CRS_TYPE *mat, const block_t *block, csr_context_t *context)
{
	SPM_CRS_IDX_TYPE nr_nzeros_blk;
	/* #non-zero elements of the block */

	SPM_CRS_IDX_TYPE nr_explicit_zeros;
	/* #explicit zeros to store */

	SPM_CRS_IDX_TYPE row_start;
	/* starting row of the block in the original matrix */

	SPM_CRS_IDX_TYPE next_row;
	SPM_CRS_IDX_TYPE next_val;

	SPM_CRS_IDX_TYPE r, c;      /* block dimensions */
	SPM_CRS_IDX_TYPE i, j, k;
	block_cache_t	*cache = context->cache;

	r = block->r;
	c = block->c;
	nr_nzeros_blk = r*c - block->pad;
	nr_explicit_zeros = block_get_nr_zeros(block);
	row_start = r*block->row;

	if (row_start + r > mat->nrows)
		/* blocked version exceeds bounds of CSR matrix */
		r = mat->nrows % r;

	if (c*block->col > mat->ncols)
		/* blocked version exceeds bounds of CSR matrix */
		c = mat->ncols % c;

	/* Skip empty rows */
	for (i = context->next_row; i <= row_start; i++)
		mat->row_ptr[i] = context->next_val;

	if (row_start < context->next_row) {
		/* The elements of the block must be inserted to existing rows */
		next_row = row_start + 1;
		next_val = mat->row_ptr[next_row];
		for (i = 0; i < r; i++) {
			SPM_CRS_IDX_TYPE nr_insert = 0;
			for (j = 0; j < c; j++) {
				SPM_CRS_IDX_TYPE nr_shift =
				    context->next_val - next_val + nr_nzeros_blk;
				SPM_CRS_IDX_TYPE where = cache->blk_map(cache->r, cache->c,
				                                        i, j);
				ELEM_TYPE value = block->values[where];
				SPM_CRS_IDX_TYPE col   = c * block->col + j;

				if (value != 0.0 || nr_explicit_zeros > 0) {
					array_insert(&mat->values[next_val], nr_shift,
					             sizeof(*mat->values), 0,
					             &value, NULL);
					array_insert(&mat->col_ind[next_val], nr_shift,
					             sizeof(*mat->col_ind), 0,
					             &col, NULL);
					next_val++;
					nr_insert++;
					if (value == 0.0)
						nr_explicit_zeros--;
				}
			}

			/* Update row_ptr */
			for (k = next_row; k <= context->next_row; k++)
				mat->row_ptr[k] += nr_insert;

			if (next_row == mat->nrows)
				/* We are at the end of the matrix */
				break;

			next_row++;
			next_val = mat->row_ptr[next_row];
		}
	} else {
		/* That's a completely new row */
		next_row = row_start;
		next_val = context->next_val;
		for (i = 0; i < r; i++, next_row++) {
			mat->row_ptr[next_row] = next_val;
			for (j = 0; j < c; j++) {
				SPM_CRS_IDX_TYPE where = cache->blk_map(cache->r, cache->c,
				                                        i, j);
				ELEM_TYPE value = block->values[where];
				SPM_CRS_IDX_TYPE col = c * block->col + j;

				if (value != 0.0 || nr_explicit_zeros > 0) {
					mat->values[next_val] = value;
					mat->col_ind[next_val] = col;
					next_val++;
					if (value == 0.0)
						nr_explicit_zeros--;
				}
			}
		}

		mat->row_ptr[next_row] = next_val;
		context->next_row = next_row;
	}

	context->next_val += nr_nzeros_blk;
	return;
}

void
csr_finalize(SPM_CRS_TYPE *mat, csr_context_t *context)
{
	SPM_CRS_IDX_TYPE i;

	for (i = context->next_row; i <= mat->nrows; i++)
		mat->row_ptr[i] = context->next_val;

	return;
}

void
bcsr_add(SPM_BCSR_TYPE *mat, const block_t *block, bcsr_context_t *context)
{
	SPM_CRS_IDX_TYPE brow, bcol, bpos;
	SPM_CRS_IDX_TYPE nr_values;
	SPM_CRS_IDX_TYPE i;

	nr_values = block->r * block->c;
	brow = block->row;
	bcol = block->col;
	bpos = context->next_block;

	for (i = context->next_brow; i <= brow; i++)
		mat->brow_ptr[i] = bpos * nr_values;

	context->next_brow = brow + 1;
	mat->bcol_ind[bpos] = bcol*mat->bc;
	memcpy(&mat->bvalues[bpos*nr_values], block->values,
	       nr_values * sizeof(ELEM_TYPE));

	context->next_block++;
	return;
}

void
bcsr_finalize(SPM_BCSR_TYPE *mat, bcsr_context_t *context)
{
	SPM_CRS_IDX_TYPE i;
	SPM_CRS_IDX_TYPE nr_values = mat->br * mat->bc;

	for (i = context->next_brow; i <= mat->nbrows; i++)
		mat->brow_ptr[i] = context->next_block * nr_values;

	return;
}
