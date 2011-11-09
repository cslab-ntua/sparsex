/*
 *  blocks.h -- Internal interface for manipulating fixed length 2-D blocks.
 *
 * Copyright (C) 2008, 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2008, 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */ 
#ifndef __BLOCKS_H__
#define __BLOCKS_H__

#include "spm_crs.h"
#include "spm_bcsr.h"
#include "bitstr.h"

enum {
    BLOCK_STORAGE_RW = 0,
    BLOCK_STORAGE_CW
};

/**
 *  A stand-alone sparse matrix element.
 */ 
typedef struct spm_elem_s {
    SPM_CRS_IDX_TYPE row;    /**< the row where the element belongs to */
    SPM_CRS_IDX_TYPE col;    /**< the column where the element belongs to */
    ELEM_TYPE value;  /**< the numerical value of the element */
} spm_elem_t;

/**
 *  A 2-D fixed length block.
 */
typedef struct block_s {
    SPM_CRS_IDX_TYPE r, c;
    /**< block dimensions */
    SPM_CRS_IDX_TYPE row;        /**< Block row */
    SPM_CRS_IDX_TYPE col;        /**< Block column */
    SPM_CRS_IDX_TYPE pad;        /**< Padding of the block */
    ELEM_TYPE *values;    /**< Block values */
    int         own;
    /**< Whether <tt>values</tt> belong to the block, i.e., have been
     * allocated during its creation. @see block_create() */
} block_t;

/**
 *  The block cache. Assumes fixed-sized <tt>r*c</tt> blocks and aligned
 *  at <tt>r</tt>- <tt>c</tt>-elements boundaries in the original matrix.
 */
typedef struct block_cache_s {
    SPM_CRS_IDX_TYPE r, c;       /**< Dimensions of blocks stored in the cache */
    SPM_CRS_IDX_TYPE size;       /**< Size of the cache (maximum #blocks allowed) */
    SPM_CRS_IDX_TYPE nr_values;
    /**< The total number of actual values stored in the cache, i.e.,
     * the size of th <tt>blk_values</tt> array. This value always equals
     * `r*c*size' */

    SPM_CRS_IDX_TYPE nr_blocks;      /**< #valid blocks */
    ELEM_TYPE *blk_values;    /**< The values of the cached blocks */
    block_t     *blocks;        /**< The cached blocks. */
    SPM_CRS_IDX_TYPE *valid_pos;
    /**< Stores the positions in the `blocks'array,
     * where valid blocks reside. */

    bitstr_t    *valid;
        /**< stores the valid bit for each position of the cache */

    SPM_CRS_IDX_TYPE (*blk_map)(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c,
                           SPM_CRS_IDX_TYPE i, SPM_CRS_IDX_TYPE j);
    /**< How 2-D blocks are mapped in memory (e.g., row-wise, column-wise) */

} block_cache_t;

/**
 *  The block cache iterator.
 */
typedef struct block_cache_iter_s {
    block_cache_t   *cache; /**< the cache to iterate on */
    SPM_CRS_IDX_TYPE     next;   /**< next element to be read */
} block_cache_iter_t;


/**
 *  The context for incrementally constructing the BCSR structure. Stores
 *  the next free position in the <tt>brow_ptr</tt> and <tt>bcol_ind</tt>
 *  arrays.
 */
typedef struct bcsr_context_s {
    SPM_CRS_IDX_TYPE next_brow;  /**< Next free position in <tt>brow_ptr</tt> */
    SPM_CRS_IDX_TYPE next_block; /**< Next free position in <tt>bcol_ind</tt> */
} bcsr_context_t;

/**
 *  The context for incrementally constructing the CSR structure. Stores
 *  the next free position in the <tt>row_ptr</tt> and <tt>values</tt>
 *  arrays.
 */
typedef struct csr_context_s {
    SPM_CRS_IDX_TYPE     next_row;   /**< Next free position in <tt>row_ptr</tt> */
    SPM_CRS_IDX_TYPE     next_val;   /**< Next free position in <tt>values</tt> */
    block_cache_t   *cache;     /**< The cache of blocks. Needed for the
                                 * retrieval of their layout */
} csr_context_t;

/**
 *  Create a new block of size <tt>r*c</tt>. Space for the values of the block
 *  will only be allocated if <tt>own</tt> is non-zero. This is done to allow
 *  for more flexible designs of elements, like block cache, that use blocks.
 *  These elements may allocate space for the values themselves and properly
 *  set up the blocks to point to the correct location.
 *
 *  @param r    The row dimension of the block.
 *  @param c    The column dimension of the block.
 *  @param own  Whether space for the values of the block should be allocated
 *              or not.
 *              
 *  @return     The new block.
 */
extern block_t *block_create(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c, int own);

/**
 *  Initializes a block. Assumes the block does not own any values and
 *  sets <tt>block->own</tt> to <tt>0</tt>.
 *
 *  @param block    The block to be initialized.
 *  @param r        The row dimension of the block.
 *  @param c        The column dimension of the block.
 */
extern void block_init(block_t *block, SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c);

/**
 *  Deletes a specific block and releases the resources held by it.
 *
 *  @param block    The block to be deleted.
 */
extern void block_delete(block_t *block);

static inline void
block_clear(block_t *block)
{
    SPM_CRS_IDX_TYPE blksize = block->r * block->c;
    memset(block->values, 0, blksize * sizeof(*(block->values)));
    block->pad = blksize;
    return;
}

/**
 *  Gets the explicit zeros stored in the block. These are NOT padding. In
 *  contrast, these elements are considered as normal elements.
 *  The reason why they exist in the block is that they also exist in the
 *  original matrix.
 *
 *  @param block    The block whose explicit zeros are queried.
 *
 *  @return         The number of explicit zeros stored in the block and are
 *                  also present in the original matrix.
 */
extern SPM_CRS_IDX_TYPE block_get_nr_zeros(const block_t *block);

/**
 *  Dump block contents to the standard error stream.
 *
 *  @param block    The block to dump.
 */
extern void block_dump(const block_t *block);

/**
 *  Creates a new block cache for <tt>size</tt> <tt>r*c</tt> blocks. The cache
 *  takes care of the proper creation and initialization of its blocks.
 *
 *  @param size The size of the new cache.
 *  @param r    The row dimension of the blocks to be stored in the cache.
 *  @param c    The column dimension of the blocks to be stored in the cache.
 *  @param storage How this block is stored in memory. This value can be
 *                 one of <tt>BLOCK_STORAGE_RW</tt> for row-wise ordering or
 *                 <tt>BLOCK_STORAGE_CW</tt> for columnt-wise ordering.
 *
 *  @return     The newly created block cache.
 */
extern block_cache_t *block_cache_create(SPM_CRS_IDX_TYPE size, SPM_CRS_IDX_TYPE r,
                                         SPM_CRS_IDX_TYPE c, int storage);

/**
 *  Deletes a block cache and releases any resources associated with it.
 *
 *  @param cache    The cache to delete.
 */
extern void block_cache_delete(block_cache_t *cache);

/**
 *  Adds a single element into the cache.
 *
 *  @param cache    The cache to insert the element into.
 *  @param elem     The element to insert into the cache.
 *
 *  @return         Whether a new block was created in the cache, to accomodate
 *  the new element.
 */
extern int block_cache_add_elem(block_cache_t *cache,
                                const spm_elem_t *elem);

/**
 *  Retrieves the total number of valid blocks in the cache.
 *
 *  @param cache    The cache whose valid blocks are queried.
 *  
 *  @return         The total number of valid blocks in the cache.
 */
static inline SPM_CRS_IDX_TYPE
block_cache_nr_blocks(block_cache_t *cache)
{
    return cache->nr_blocks;
}

/**
 *  Clears the cache. It simply invalidates and zeros out all the valid
 *  blocks inside the cache.
 *
 *  @param  cache   The cache to clear.
 */
extern void block_cache_clear(block_cache_t *cache);


/**
 *  Sorts the valid block of cache.
 *
 *  @param cache    The cache.
 */ 
extern void block_cache_sort(block_cache_t *cache);

/**
 *  Creates a new iterator for a cache.
 *
 *  @param cache    The cache for which an iterator will be created.
 *  
 *  @return         The new iterator.
 */
extern block_cache_iter_t *block_cache_iter_create(block_cache_t *cache);

/**
 *  Attaches an iterator to a specific cache.
 *
 *  @param iter     The iterator to attach.
 *  @param cache    The cache to attach the iterator to.
 */
static inline void
block_cache_iter_init(block_cache_iter_t *iter, block_cache_t *cache)
{
    iter->cache = cache;
    iter->next  = 0;
    return;
}

/**
 *  Positions the iterator to the beginning of the cache.
 *
 *  @param iter     The iterator to position.
 */ 
static inline void
block_cache_iter_begin(block_cache_iter_t *iter)
{
    iter->next = 0;
    return;
}

/**
 *  Positions the iterator to the end of the cache.
 *
 *  @param iter     The iterator to position.
 */
static inline void
block_cache_iter_end(block_cache_iter_t *iter)
{
    iter->next = iter->cache->nr_blocks;
    return;
}

/**
 *  Checks whether more valid blocks exist in the cache to iterate on
 *
 *  @param iter     The iterator.
 */ 
static inline int
block_cache_iter_more(block_cache_iter_t *iter)
{
    return (iter->next < iter->cache->nr_blocks);
}

/**
 *  Returns the next valid block in the cache and proceeds with the iteration.
 *  
 *  @param iter     The iterator.
 *
 *  @return         A pointer to the next block in the iteration.
 */
static inline void
block_cache_iter_next(block_cache_iter_t *iter, block_t *block)
{
    block_cache_t   *cache = iter->cache;

    *block = cache->blocks[cache->valid_pos[iter->next++]];
    return;
}

/**
 *  Adds the elements of `block' to `mat' in CSR format.
 */
extern void csr_add(SPM_CRS_TYPE *mat, const block_t *block,
                    csr_context_t *context);

extern void bcsr_add(SPM_BCSR_TYPE *mat, const block_t *block,
                     bcsr_context_t *context);

extern void bcsr_finalize(SPM_BCSR_TYPE *mat, bcsr_context_t *context);

extern void csr_finalize(SPM_CRS_TYPE *mat, csr_context_t *context);

/**
 *  Storage layouts of blocks assumed by specific block multiplication
 *  routines.
 */
struct storage_layout_s {
    SPM_CRS_IDX_TYPE sl_r;       /**< #rows of the block */
    SPM_CRS_IDX_TYPE sl_c;       /**< #columns of the block */
    int         sl_simd;    /**< is an SIMD impl */
    int         sl_layout;  /**< the assumed storage layout */
};

extern struct storage_layout_s blk_layout_list[];
extern const int blk_layout_default;
extern const int blk_layout_list_size;
extern const int blk_layout_list_begin;
extern const int blk_layout_list_end;

#endif  /* __BLOCKS_H__ */
