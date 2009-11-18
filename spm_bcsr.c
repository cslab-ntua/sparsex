#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "spm_crs.h"
#include "spm_mt.h"
#include "spm_bcsr.h"
#include "mt_lib.h"

#ifndef SPM_CRS_BITS
#define SPM_CRS_BITS 32
#endif

#ifndef ELEM_TYPE
#define ELEM_TYPE double
#endif

#define BCSR_DEFAULT_COL 2
#define BCSR_DEFAULT_ROW 2

///////////////////////////////////////////////////////////////////////////////
// block.{h,c} START
///////////////////////////////////////////////////////////////////////////////

/* Ceiling of two integers */
#define iceil(a, b) ((a) / (b) + ((a) % (b) != 0))

#include <limits.h>

enum {
	BLOCK_STORAGE_RW = 0,
	BLOCK_STORAGE_CW
};


void *
Malloc(size_t num)
{
    void *new = malloc(num);
    if (!new){
		perror("malloc");
        exit(1);
	}
    return new;
}

void *
Calloc(size_t num, size_t size)
{
    void *new = calloc(num, size);
	if (!new){
		perror("malloc");
		exit(1);
	}
    return new;
}

#define SPM_CALLOC(type, num)                               \
        ((type *) Calloc ((num), sizeof(type)))
#define SPM_MALLOC(type, num)                               \
        ((type *) Malloc ((num) * sizeof(type)))
#define SPM_FREE(stale)                         do {        \
        if (stale) { free (stale);  stale = 0; }            \
                                                } while (0)

#ifndef LONG_BIT
#   define LONG_BIT_DEFINED 1
#   define LONG_BIT (CHAR_BIT*sizeof(long))
#endif

typedef struct bitstr_s {
    unsigned long   *base_addr; /* the base address of the bit-string */
    unsigned long   nr_bits;    /* the size of the bit-string */
} bitstr_t;

/* Get the actual word in the bitstr to compare with */
#define BS_WORD(pos, word, word_pos, word_cmp) { \
        *word     = (unsigned long) pos / LONG_BIT;  \
        *word_pos = (unsigned long) pos % LONG_BIT;  \
        *word_cmp = 1UL << *word_pos;                \
    }

inline static bitstr_t *
bs_create(unsigned long nr_bits)
{
    bitstr_t        *new = SPM_MALLOC(bitstr_t, 1);
    unsigned long   nr_words = iceil(nr_bits, LONG_BIT);

    new->base_addr = SPM_CALLOC(unsigned long, nr_words);
    new->nr_bits   = nr_words * LONG_BIT;
    return new;
}

inline static void
bs_delete(bitstr_t *bs)
{
    if (!bs)
        return;

    SPM_FREE(bs->base_addr);
    SPM_FREE(bs);
    return;
}

inline static int
bs_isset(const bitstr_t *bs, unsigned long pos)
{
    unsigned long   word;
    unsigned long   word_cmp;
    unsigned long   word_pos;

    BS_WORD(pos, &word, &word_pos, &word_cmp);
    return ((bs->base_addr[word] & word_cmp) == word_cmp);
}

inline static void
bs_set(bitstr_t *bs, unsigned long pos)
{
    unsigned long   word;
    unsigned long   word_cmp;
    unsigned long   word_pos;

    BS_WORD(pos, &word, &word_pos, &word_cmp);
    bs->base_addr[word] |= word_cmp;
//    printf("bs_set:base_addr = 0x%lx\n", (unsigned long) bs->base_addr);
    return;
}

inline static void
bs_clear(bitstr_t *bs, unsigned long pos)
{
    unsigned long   word;
    unsigned long   word_cmp;
    unsigned long   word_pos;

    BS_WORD(pos, &word, &word_pos, &word_cmp);
    bs->base_addr[word] &= ~word_cmp;
    return;
}

inline static void
bs_setall(bitstr_t *bs)
{
    memset(bs->base_addr, 0xff, iceil(bs->nr_bits, CHAR_BIT));
    return;
}

inline static void
bs_clearall(bitstr_t *bs)
{
//    printf("bs_clearall:in:base_addr = 0x%lx\n", (unsigned long) bs->base_addr);
    memset(bs->base_addr, 0, iceil(bs->nr_bits, CHAR_BIT));
    return;
}

inline static unsigned long
bs_size(bitstr_t *bs)
{
    return bs->nr_bits;
}

inline static void
bs_print(bitstr_t *bs)
{
    unsigned long   nr_words = bs->nr_bits / LONG_BIT;

    unsigned long   i;
    printf("0x");
    /* Print most significant word first */
    for (i = nr_words - 1; i != 0; i--)
        printf("%016lx ", bs->base_addr[i]);

    /* Print the `i == 0' case explicitly */
    printf("%016lx", bs->base_addr[0]);
    printf("\n");
}

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

static inline void
block_clear(block_t *block)
{
    SPM_CRS_IDX_TYPE blksize = block->r * block->c;
    memset(block->values, 0, blksize * sizeof(*(block->values)));
    block->pad = blksize;
    return;
}

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

#if 0
/**
 *  Storage layouts of blocks assumed by specific block multiplication
 *  routines.
 */
static struct storage_layout_s {
    SPM_CRS_IDX_TYPE sl_r;       /**< #rows of the block */
    SPM_CRS_IDX_TYPE sl_c;       /**< #columns of the block */
    int              sl_simd;    /**< is an SIMD impl */
    int              sl_layout;  /**< the assumed storage layout */
} blk_layout_list[] = {
    {
        .sl_r = 0,
        .sl_c = 0,
        .sl_layout = BLOCK_STORAGE_RW,
        .sl_simd = 0
    },    /* The default layout */

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
        .sl_layout = ((SPM_VALUE_BINDING == 'd') ?
                      BLOCK_STORAGE_CW : BLOCK_STORAGE_RW),
        .sl_simd = 1
    },

    {
        .sl_r = 4,
        .sl_c = 2,
        .sl_layout = BLOCK_STORAGE_CW,
        .sl_simd = 1
    }
}

#define BLK_LAYOUT_LIST_SIZE                           \
    (sizeof(blk_layout_list) / sizeof(*blk_layout_list))

const int blk_layout_default = 0;
const int blk_layout_list_size = BLK_LAYOUT_LIST_SIZE;
const int blk_layout_list_begin = 0;
const int blk_layout_list_end = BLK_LAYOUT_LIST_SIZE - 1;
#endif

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

static SPM_CRS_IDX_TYPE
map_cw(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c, SPM_CRS_IDX_TYPE i, SPM_CRS_IDX_TYPE j)
{
    return j*r + i;
}

static SPM_CRS_IDX_TYPE
map_rw(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c, SPM_CRS_IDX_TYPE i, SPM_CRS_IDX_TYPE j)
{
    return i*c + j;
}

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
block_t *
block_create(SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c, int own)
{
    block_t  *new = SPM_MALLOC(block_t, 1);

    new->r   = r;
    new->c   = c;
    //new->row = SPM_INDEX_INVALID;
    //new->col = SPM_INDEX_INVALID;
    new->pad = r*c; /* the block is full of zeros */
    new->own = own;

    if (own)
        /* Allocate your own values. */
        new->values = SPM_CALLOC(ELEM_TYPE, r*c);
    else
        new->values = NULL;

    return new;
}

/**
 *  Initializes a block. Assumes the block does not own any values and
 *  sets <tt>block->own</tt> to <tt>0</tt>.
 *
 *  @param block    The block to be initialized.
 *  @param r        The row dimension of the block.
 *  @param c        The column dimension of the block.
 */
void
block_init(block_t *block, SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c)
{
    block->r = r;
    block->c = c;
    //block->row = SPM_INDEX_INVALID;
    //block->col = SPM_INDEX_INVALID;
    block->pad = r*c;
    block->own = 0;
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
uint64_t
block_get_nr_zeros(const block_t *block)
{
    uint64_t nr_zeros;
    uint64_t i;

    for (i = 0, nr_zeros = 0; i < block->r * block->c; i++)
        if (block->values[i] == 0.0)
            nr_zeros++;

    return (nr_zeros - block->pad);
}

/**
 *  Deletes a specific block and releases the resources held by it.
 *
 *  @param block    The block to be deleted.
 */
void
block_delete(block_t *block)
{
    if (block && block->own)
        SPM_FREE(block->values);

    SPM_FREE(block);
    return;
}

/**
 *  Dump block contents to the standard error stream.
 *
 *  @param block    The block to dump.
 *
 *  FIXME: Should use the `SPM_INDEX_FMT' and `SPM_VALUE_FMT' constants
 *  for the output format.
 */
#if 0
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
 *  Block cache operations.
 */

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
static block_cache_t *
block_cache_create(SPM_CRS_IDX_TYPE size, SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c, int storage)
{
    block_cache_t   *new;

    new = SPM_MALLOC(block_cache_t, 1);
    new->blocks     = SPM_MALLOC(block_t, size);
    /* values are initialized to zero */
    new->blk_values = SPM_CALLOC(ELEM_TYPE, size*r*c);
    new->valid      = bs_create(size);
    new->valid_pos  = SPM_CALLOC(SPM_CRS_IDX_TYPE, size);
    new->r          = r;
    new->c          = c;
    new->size       = size;
    new->nr_values  = r*c*size;
    new->nr_blocks  = 0;

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

/**
 *  Deletes a block cache and releases any resources associated with it.
 *
 *  @param cache    The cache to delete.
 */
static void
block_cache_delete(block_cache_t *cache)
{
    if (!cache)
        return;

    bs_delete(cache->valid);
    SPM_FREE(cache->valid_pos);
    SPM_FREE(cache->blocks);
    SPM_FREE(cache->blk_values);
    SPM_FREE(cache);
    return;
}

/**
 *  Sorts the valid block of cache.
 *
 *  @param cache    The cache.
 */
static void
block_cache_sort(block_cache_t *cache)
{
    qsort(cache->valid_pos, cache->nr_blocks, sizeof(SPM_CRS_IDX_TYPE),
          spm_index_cmp);

    return;
}

/**
 *  Clears the cache. It simply invalidates and zeros out all the valid
 *  blocks inside the cache.
 *
 *  @param  cache   The cache to clear.
 */
static void
block_cache_clear(block_cache_t *cache)
{
    SPM_CRS_IDX_TYPE i;

    for (i = 0; i < cache->nr_blocks; i++)
        bs_clear(cache->valid, cache->valid_pos[i]);

    cache->nr_blocks = 0;
    return;
}

/**
 *  Adds a single element into the cache.
 *
 *  @param cache    The cache to insert the element into.
 *  @param elem     The element to insert into the cache.
 *
 *  @return         Whether a new block was created in the cache, to accomodate
 *  the new element.
 */
static int
block_cache_add_elem(block_cache_t *cache, const spm_elem_t *elem)
{
    SPM_CRS_IDX_TYPE bcol;
    SPM_CRS_IDX_TYPE brow;
    SPM_CRS_IDX_TYPE br, bc;
    int         is_new;

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

    if (cache->blocks[bcol].pad <= 0) {
        //DEBUG_PRINT_LONG(bcol);
        //block_dump(&cache->blocks[bcol]);
    }

    assert(cache->blocks[bcol].pad);
    if (cache->blocks[bcol].pad != 0)
        cache->blocks[bcol].pad--;

    return is_new;
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

    for (i = context->next_brow; i <= mat->nr_brows; i++)
        mat->brow_ptr[i] = context->next_block * nr_values;

    return;
}

///////////////////////////////////////////////////////////////////////////////
// block.{c,h} END
///////////////////////////////////////////////////////////////////////////////


static void
_scan_csr(SPM_CRS_TYPE *crs, int r, int c,
          block_cache_t *blk_cache, uint64_t *nr_blocks)
{
    SPM_CRS_IDX_TYPE i, j, k;
    spm_elem_t  curr_elem;

    *nr_blocks    = 0;
    //*nr_leftovers = 0;
    for (i = 0; i < crs->nrows; i += r) {
        for (j = 0; j < r; j++) {
            curr_elem.row = i + j;
            if (curr_elem.row + 1 > crs->nrows)
                break;
            for (k = crs->row_ptr[i + j]; k < crs->row_ptr[i + j + 1]; k++) {
                curr_elem.col   = crs->col_ind[k];
                curr_elem.value = crs->values[k];
                block_cache_add_elem(blk_cache, &curr_elem);
            }
        }

        *nr_blocks += block_cache_nr_blocks(blk_cache);

		#if 0
        block_cache_iter_t  iter;
        block_t             block;
        /* Skip the search for leftovers */
		skip = 1;
		long max_pad;
        if (!skip) {
            block_cache_iter_init(&iter, blk_cache);
            block_cache_iter_begin(&iter);
            while (block_cache_iter_more(&iter)) {
                block_cache_iter_next(&iter, &block);
                if (block.pad > max_pad) {
                    *nr_leftovers += r*c - block.pad;
                    (*nr_blocks)--;
                }
            }
        }
		#endif

        block_cache_clear(blk_cache);   /* clear the cache */
    }

   return;
}

static void
_fill_bcsr(SPM_BCSR_TYPE *bcsr, SPM_CRS_TYPE *crs, block_cache_t *blk_cache)
{
    SPM_CRS_IDX_TYPE     i, j, k;
    spm_elem_t      curr_elem;
    SPM_CRS_IDX_TYPE     r;
    bcsr_context_t  bcsr_context;

    bcsr_context.next_brow  = 0;
    bcsr_context.next_block = 0;

	#if 0
    csr_context_t   csr_context;
    csr_context.next_row    = 0;
    csr_context.next_val    = 0;
    csr_context.cache       = blk_cache;
	#endif

    r = bcsr->br;
    for (i = 0; i < crs->nrows; i += r) {
        for (j = 0; j < r; j++) {
            curr_elem.row = i + j;
            if (curr_elem.row + 1 > crs->nrows)
                break;
            for (k = crs->row_ptr[i + j];
                 k < crs->row_ptr[i + j + 1]; k++) {
                curr_elem.col   = crs->col_ind[k];
                curr_elem.value = crs->values[k];
                block_cache_add_elem(blk_cache, &curr_elem);
            }
        }

        block_cache_iter_t  iter;
        block_t             block;

        /*
         * First sort the valid blocks. This is not a necessary step; just to
         * be strictly conformant to the BCSR format.
         */
        block_cache_sort(blk_cache);
        block_cache_iter_init(&iter, blk_cache);
        block_cache_iter_begin(&iter);
        while (block_cache_iter_more(&iter)) {
            block_cache_iter_next(&iter, &block);
			#if 0
            if (block.pad > max_pad) {
                /* This block is added to the leftover */
                if (leftover)
                    csr_add(leftover, &block, &csr_context);
            } else {
                bcsr_add(mat, &block, &bcsr_context);
            }
			#endif
			bcsr_add(bcsr, &block, &bcsr_context);
        }

        block_cache_clear(blk_cache);   /* clear the cache */
    }

    bcsr_finalize(bcsr, &bcsr_context);
	#if 0
    if (leftover)
        csr_finalize(leftover, &csr_context);
	#endif

    return;
}

static SPM_BCSR_TYPE *
crs_to_bcsr(SPM_CRS_TYPE *crs, int r, int c)
{
	SPM_BCSR_TYPE *bcsr;
	block_cache_t *blk_cache;
	uint64_t blk_cache_size, nr_blocks;

	bcsr = malloc(sizeof(SPM_BCSR_TYPE));
	if (!bcsr){
		perror("malloc");
		exit(1);
	}

	#if 0
	char *tmp;
	int c, r;
	r = (tmp = getenv("BCSR_ROW")) ? atol(tmp) : BCSR_DEFAULT_ROW;
	c = (tmp = getenv("BCSR_COL")) ? atol(tmp) : BCSR_DEFAULT_COL;
	assert (c != 0 && r != 0);
	printf("bcsr r:%ld c:%ld\n", r, c);
	#endif

	// scan
	blk_cache_size = iceil(crs->ncols, c);
	blk_cache = block_cache_create(blk_cache_size, r, c, BLOCK_STORAGE_RW);
	_scan_csr(crs, r, c, blk_cache, &nr_blocks);
	block_cache_delete(blk_cache);

	// create
    bcsr->nr_brows  = iceil(crs->nrows, r);
    bcsr->brow_ptr  = SPM_MALLOC(SPM_CRS_IDX_TYPE, bcsr->nr_brows + 1);
    bcsr->bcol_ind  = SPM_MALLOC(SPM_CRS_IDX_TYPE, nr_blocks);
    bcsr->bvalues   = SPM_MALLOC(ELEM_TYPE, nr_blocks*r*c);
    /* If got so far, everything is ok; continue with init'ion and exit. */
    bcsr->nr_nzeros = crs->nz;
    bcsr->nr_blocks = nr_blocks;
    bcsr->nr_rows   = bcsr->nr_brows * r;
    bcsr->nr_cols   = iceil(crs->ncols, c) * c;
    bcsr->br        = r;
    bcsr->bc        = c;
    bcsr->storage   = BLOCK_STORAGE_RW;
	#if 0
    bcsr->storage   = blk_layout_list[blk_layout_default].sl_layout;
    /* Search for a special storage layout to be used. */
    int i;
    for (i = blk_layout_list_begin; i < blk_layout_list_size; i++) {
        if (blk_layout_list[i].sl_simd == SPM_USE_SIMD  &&
            blk_layout_list[i].sl_r == r                &&
            blk_layout_list[i].sl_c == c) {
            bcsr->storage = blk_layout_list[i].sl_layout;
            break;
        }
    }
	#endif

	// init
	blk_cache_size = iceil(bcsr->nr_cols, bcsr->bc);
	blk_cache = block_cache_create(blk_cache_size, bcsr->br, bcsr->bc, bcsr->storage);
	_fill_bcsr(bcsr, crs, blk_cache);
	block_cache_delete(blk_cache);

	return bcsr;
}

static void split_rows(SPM_BCSR_TYPE *mat,
                       SPM_CRS_IDX_TYPE nr_splits,
                       SPM_CRS_IDX_TYPE *splits)
{

    SPM_CRS_IDX_TYPE nr_nzeros = mat->nr_blocks * mat->br * mat->bc;
    SPM_CRS_IDX_TYPE *brow_ptr = mat->brow_ptr;
    SPM_CRS_IDX_TYPE min_nzeros_per_split = nr_nzeros / nr_splits;
    SPM_CRS_IDX_TYPE i;

    SPM_CRS_IDX_TYPE nr_nzeros_per_split = 0;
    SPM_CRS_IDX_TYPE split_pos = 1;


    splits[0] = 0;
    /* Iterate over block rows */
    for (i = 0; i < mat->nr_brows; i++) {
        nr_nzeros_per_split += brow_ptr[i+1] - brow_ptr[i];
        if (nr_nzeros_per_split > min_nzeros_per_split) {
            /* new split point found */
            splits[split_pos++] = (i+1)*mat->br;
            nr_nzeros_per_split -= min_nzeros_per_split;
        }
    }

    splits[nr_splits] = mat->nr_brows * mat->br;
    return;
}


static void *
init_mmf_wrap(char *mmf_file,
              unsigned long *rows_nr, unsigned long *cols_nr,
              unsigned long *nz_nr, int r, int c)
{
	spm_mt_t *spm_mt;
	SPM_CRS_TYPE *crs;
	SPM_BCSR_TYPE *bcsr;
	SPM_BCSR_MT_TYPE *bcsr_mts;
	SPM_CRS_IDX_TYPE *splits;
	unsigned int nr_cpus, *cpus;
	long i;

	mt_get_options(&nr_cpus, &cpus);
	spm_mt = malloc(sizeof(spm_mt_t));
	if (!spm_mt){
		perror("malloc");
		exit(1);
	}

	spm_mt->nr_threads = nr_cpus;
	spm_mt->spm_threads = malloc(sizeof(spm_mt_thread_t)*nr_cpus);
	bcsr_mts = malloc(sizeof(SPM_BCSR_MT_TYPE)*nr_cpus);
	splits = malloc(sizeof(SPM_CRS_IDX_TYPE)*(nr_cpus + 1));
	if (!spm_mt->spm_threads || !bcsr_mts || !splits){
		perror("malloc");
		exit(1);
	}
	crs = SPM_CRS_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);
	bcsr = crs_to_bcsr(crs, r, c);
	split_rows(bcsr, nr_cpus, splits);

	for (i=0; i < spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread;
		SPM_BCSR_MT_TYPE *bcsr_mt;

		spm_thread = spm_mt->spm_threads + i;
		bcsr_mt = bcsr_mts + i;

		bcsr_mt->row_start = splits[i];
		bcsr_mt->row_end = splits[i+1];
		bcsr_mt->nnz_nr = bcsr->brow_ptr[bcsr_mt->row_end/r] - bcsr->brow_ptr[bcsr_mt->row_start/r];
		bcsr_mt->bcsr = bcsr;
		//printf("bcsr_mt: thread:%ld row_start:%lu row_end:%lu\n", i, bcsr_mt->row_start, bcsr_mt->row_end);

		spm_thread->spm = bcsr_mt;
		spm_thread->spmv_fn = NULL;
	}

	free(splits);
	SPM_CRS_NAME(_destroy)(crs);
	free(cpus);
	return spm_mt;
}


static void inline bcsr_spmv_wrap(void *spm, VECTOR_TYPE * __restrict__ in, VECTOR_TYPE *__restrict__ out, int r, int c)
{
	SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
	SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
	SPM_CRS_IDX_TYPE i, _i, j, _j, k, l;

     	ELEM_TYPE * __restrict__ bvalues = bcsr->bvalues;
     	ELEM_TYPE * __restrict__ x = in->elements;
     	ELEM_TYPE * __restrict__ y = out->elements;
	SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
	SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
	SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
	SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

	for (i = r_start, _i = r_start / r; i < r_end; i += r, _i++) {
		for (j = brow_ptr[_i], _j = j / (r*c); j < brow_ptr[_i+1];
			 j += r*c, _j++) {
			SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
			for (k = 0; k < r; k++) {
				//y[i+k] = 0;
				for (l = 0; l < c; l++){
					SPM_CRS_IDX_TYPE y_indx;
					SPM_CRS_IDX_TYPE x_indx;
					SPM_CRS_IDX_TYPE v_indx;
					x_indx = x_start + l;
					y_indx = i + k;
					v_indx = j+k*c + l;
					//y[i+k] += bvalues[j+k*c+l]*x[x_start + l];
					//printf("y[%lu] += bvalues[%lu] * x[%lu]\n", y_indx, v_indx, x_indx);
					//printf("%lf += %lf * %lf \n", y[y_indx], bvalues[v_indx], x[x_indx]);
					y[y_indx] += bvalues[v_indx] * x[x_indx];
				}
			}
		}
	}

	return;
}

static uint64_t bcsr_size(void *spm)
{
	spm_mt_t *spm_mt;
	SPM_BCSR_TYPE *bcsr;
	unsigned long r, c, ret = 0;

	spm_mt = spm;
	bcsr = ((SPM_BCSR_MT_TYPE *)spm_mt->spm_threads->spm)->bcsr;
	r = bcsr->br;
	c = bcsr->bc;
	ret +=  bcsr->nr_blocks*r*c*sizeof(ELEM_TYPE);
	ret += (bcsr->nr_blocks)*sizeof(SPM_CRS_IDX_TYPE);
	ret += (bcsr->nr_brows + 1)*sizeof(SPM_CRS_IDX_TYPE);
	return ret;
}

static void bcsr_destroy(void *spm)
{
	return;
}

#define DECLARE_BCSR(c, r)  \
static void *\
bcsr_init_mmf_ ## r ## x ## c(char *f,                   \
                              unsigned long *rows,       \
                              unsigned long *cols,       \
                              unsigned long *nnz)        \
{                                                        \
	return init_mmf_wrap(f, rows, cols, nnz, r, c);  \
}                                                        \
                                                         \
static void spmv_bcsr_ ## r ## x ## c (void *spm,        \
                                       VECTOR_TYPE *in,  \
                                       VECTOR_TYPE *out) \
{                                                        \
       return bcsr_spmv_wrap(spm, in, out, r, c);        \
}                                                        \
XSPMV_MT_METH_INIT(spmv_bcsr_ ## r ## x ## c,            \
               bcsr_init_mmf_ ## r ## x ## c,            \
               bcsr_size,                                \
               bcsr_destroy,                             \
               sizeof(ELEM_TYPE))                        \

DECLARE_BCSR(1,1)
DECLARE_BCSR(1,2)
DECLARE_BCSR(1,3)
DECLARE_BCSR(1,4)
DECLARE_BCSR(2,1)
DECLARE_BCSR(2,2)
DECLARE_BCSR(2,3)
DECLARE_BCSR(2,4)
DECLARE_BCSR(3,1)
DECLARE_BCSR(3,2)
DECLARE_BCSR(3,3)
DECLARE_BCSR(3,4)
DECLARE_BCSR(4,1)
DECLARE_BCSR(4,2)
DECLARE_BCSR(4,3)
DECLARE_BCSR(4,4)

void
bcsr_mulv_1x2(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, j, _j;

    for (i = r_start; i < r_end; i++) {
        register ELEM_TYPE    _y0 = 0;
        for (j = brow_ptr[i], _j = j / (1*2); j < brow_ptr[i+1];
             j += 1*2, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j] * x[x_start] + bvalues[j+1] * x[x_start + 1];
        }

        y[i] = _y0;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_1x2, bcsr_init_mmf_1x2, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))



void
bcsr_mulv_1x3(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, j, _j;

    for (i = r_start; i < r_end; i++) {
        register ELEM_TYPE    _y0 = 0;
        for (j = brow_ptr[i], _j = j / (1*3); j < brow_ptr[i+1];
             j += 1*3, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j] * x[x_start] +
                bvalues[j+1] * x[x_start + 1] +
                bvalues[j+2] * x[x_start + 2];
        }

        y[i] = _y0;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_1x3, bcsr_init_mmf_1x3, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))



void
bcsr_mulv_1x4(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, j, _j;

    for (i = r_start; i < r_end; i++) {
        register ELEM_TYPE    _y0 = 0;
        for (j = brow_ptr[i], _j = j / (1*4); j < brow_ptr[i+1];
             j += 1*4, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j] * x[x_start] +
                bvalues[j+1] * x[x_start + 1] +
                bvalues[j+2] * x[x_start + 2] +
                bvalues[j+3] * x[x_start + 3];
        }

        y[i] = _y0;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_1x4, bcsr_init_mmf_1x4, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))


void
bcsr_mulv_2x1(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 2; i < r_end; i += 2, _i++) {
        register ELEM_TYPE _y0 = 0;
        register ELEM_TYPE _y1 = 0;
        for (j = brow_ptr[_i], _j = j / (2*1); j < brow_ptr[_i+1];
             j += 2*1, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start];
            _y1 += bvalues[j+1] * x[x_start];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_2x1, bcsr_init_mmf_2x1, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))


void
bcsr_mulv_2x2(void *spm, VECTOR_TYPE * __restrict__ in, VECTOR_TYPE * __restrict__ out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;
     SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 2; i < r_end; i += 2, _i++) {
        register ELEM_TYPE _y0 = 0;
        register ELEM_TYPE _y1 = 0;
        for (j = brow_ptr[_i], _j = j / (2*2); j < brow_ptr[_i+1];
             j += 2*2, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start] + bvalues[j+1] * x[x_start+1];
            _y1 += bvalues[j+2] * x[x_start] + bvalues[j+3] * x[x_start+1];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
    }

    return;
}
XSPMV_MT_METH_INIT(bcsr_mulv_2x2, bcsr_init_mmf_2x2, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))

void
bcsr_mulv_2x3(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 2; i < r_end; i += 2, _i++) {
        register ELEM_TYPE _y0 = 0;
        register ELEM_TYPE _y1 = 0;
        for (j = brow_ptr[_i], _j = j / (2*3); j < brow_ptr[_i+1];
             j += 2*3, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start] + bvalues[j+1] * x[x_start+1] +
                bvalues[j+2] * x[x_start+2];
            _y1 += bvalues[j+3] * x[x_start] + bvalues[j+4] * x[x_start+1] +
                bvalues[j+5] * x[x_start+2];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_2x3, bcsr_init_mmf_2x3, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))



void
bcsr_mulv_2x4(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 2; i < r_end; i += 2, _i++) {
        register ELEM_TYPE _y0 = 0;
        register ELEM_TYPE _y1 = 0;
        for (j = brow_ptr[_i], _j = j / (2*4); j < brow_ptr[_i+1];
             j += 2*4, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start]   + bvalues[j+1] * x[x_start+1] +
                bvalues[j+2] * x[x_start+2] + bvalues[j+3] * x[x_start+3];
            _y1 += bvalues[j+4] * x[x_start]   + bvalues[j+5] * x[x_start+1] +
                bvalues[j+6] * x[x_start+2] + bvalues[j+7] * x[x_start+3];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_2x4, bcsr_init_mmf_2x4, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))



void
bcsr_mulv_3x1(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 3; i < r_end; i += 3, _i++) {
        register ELEM_TYPE _y0 = 0;
        register ELEM_TYPE _y1 = 0;
        register ELEM_TYPE _y2 = 0;
        for (j = brow_ptr[_i], _j = j / (3*1); j < brow_ptr[_i+1];
             j += 3*1, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start];
            _y1 += bvalues[j+1] * x[x_start];
            _y2 += bvalues[j+2] * x[x_start];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
        y[i+2] = _y2;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_3x1, bcsr_init_mmf_3x1, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))



void
bcsr_mulv_3x2(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 3; i < r_end; i += 3, _i++) {
        register ELEM_TYPE _y0 = 0;
        register ELEM_TYPE _y1 = 0;
        register ELEM_TYPE _y2 = 0;
        for (j = brow_ptr[_i], _j = j / (3*2); j < brow_ptr[_i+1];
             j += 3*2, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start] + bvalues[j+1] * x[x_start+1];
            _y1 += bvalues[j+2] * x[x_start] + bvalues[j+3] * x[x_start+1];
            _y2 += bvalues[j+4] * x[x_start] + bvalues[j+5] * x[x_start+1];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
        y[i+2] = _y2;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_3x2, bcsr_init_mmf_3x2, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))



void
bcsr_mulv_4x1(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 4; i < r_end; i += 4, _i++) {
        ELEM_TYPE _y0 = 0;
        ELEM_TYPE _y1 = 0;
        ELEM_TYPE _y2 = 0;
        ELEM_TYPE _y3 = 0;
        for (j = brow_ptr[_i], _j = j / (4*1); j < brow_ptr[_i+1];
             j += 4*1, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start];
            _y1 += bvalues[j+1] * x[x_start];
            _y2 += bvalues[j+2] * x[x_start];
            _y3 += bvalues[j+3] * x[x_start];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
        y[i+2] = _y2;
        y[i+3] = _y3;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_4x1, bcsr_init_mmf_4x1, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))



void
bcsr_mulv_4x2(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
     SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm;
     SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
     ELEM_TYPE *bvalues = bcsr->bvalues;
     ELEM_TYPE *x = in->elements;
     ELEM_TYPE *y = out->elements;
     SPM_CRS_IDX_TYPE *brow_ptr = bcsr->brow_ptr;
     SPM_CRS_IDX_TYPE *bcol_ind = bcsr->bcol_ind;
     SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
     SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

    SPM_CRS_IDX_TYPE i, _i, j, _j;

    for (i = r_start, _i = r_start / 4; i < r_end; i += 4, _i++) {
        ELEM_TYPE _y0 = 0;
        ELEM_TYPE _y1 = 0;
        ELEM_TYPE _y2 = 0;
        ELEM_TYPE _y3 = 0;
        for (j = brow_ptr[_i], _j = j / (4*2); j < brow_ptr[_i+1];
             j += 4*2, _j++) {
            SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
            _y0 += bvalues[j]   * x[x_start] + bvalues[j+1] * x[x_start+1];
            _y1 += bvalues[j+2] * x[x_start] + bvalues[j+3] * x[x_start+1];
            _y2 += bvalues[j+4] * x[x_start] + bvalues[j+5] * x[x_start+1];
            _y3 += bvalues[j+6] * x[x_start] + bvalues[j+7] * x[x_start+1];
        }

        y[i]   = _y0;
        y[i+1] = _y1;
        y[i+2] = _y2;
        y[i+3] = _y3;
    }

    return;

}
XSPMV_MT_METH_INIT(bcsr_mulv_4x2, bcsr_init_mmf_4x2, bcsr_size, bcsr_destroy, sizeof(ELEM_TYPE))

#ifdef SPM_NUMA

//#include <numa.h>

static int numa_node_from_cpu(int cpu)
{
#if 0
	struct bitmask *bmp;
	int nnodes, node, ret;

	bmp = numa_allocate_cpumask();
	nnodes =  numa_num_configured_nodes();
	for (node = 0; node < nnodes; node++){
		numa_node_to_cpus(node, bmp);
		if (numa_bitmask_isbitset(bmp, cpu)){
			ret = node;
			goto end;
		}
	}
	ret = -1;
end:
	numa_bitmask_free(bmp);
	return ret;
#endif
	return 0;
}

void *numa_alloc_onnode(unsigned long size, int node)
{
	return malloc(size);
}

static void *
init_mmf_numa_wrap(char *mmf_file,
              unsigned long *rows_nr, unsigned long *cols_nr,
              unsigned long *nz_nr, int r, int c)
{
	spm_mt_t *spm_mt;
	SPM_BCSR_TYPE *bcsr;
	SPM_CRS_IDX_TYPE *brow_ptr;
	SPM_CRS_IDX_TYPE *bcol_ind;
	ELEM_TYPE *bvalues;

	spm_mt = init_mmf_wrap(mmf_file, rows_nr, cols_nr, nz_nr, r, c);

	bcsr = ((SPM_BCSR_MT_TYPE *)spm_mt->spm_threads->spm)->bcsr;
	brow_ptr = bcsr->brow_ptr;
	bcol_ind = bcsr->bcol_ind;
	bvalues = bcsr->bvalues;

	unsigned long i, total_nnz=0;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm_thread->spm;
		/* get numa node from cpu */
		int node = numa_node_from_cpu(spm_thread->cpu);
		/* allocate new space */
		SPM_BCSR_TYPE *numa_bcsr = numa_alloc_onnode(sizeof(SPM_BCSR_TYPE), node);
		if (!numa_bcsr){
			perror("numa_alloc_onnode");
			exit(1);
		}
		uint64_t nnz = bcsr_mt->nnz_nr;
		uint64_t nrows = bcsr_mt->row_end - bcsr_mt->row_start;
		numa_bcsr->bvalues = numa_alloc_onnode(sizeof(ELEM_TYPE)*nnz, node);
		numa_bcsr->bcol_ind = numa_alloc_onnode(sizeof(SPM_CRS_IDX_TYPE)*nnz/(r*c), node);
		numa_bcsr->brow_ptr = numa_alloc_onnode(sizeof(SPM_CRS_IDX_TYPE)*(nrows/r+1), node);
		if (!numa_bcsr->bvalues || !numa_bcsr->bcol_ind | !numa_bcsr->brow_ptr){
			perror("numa_alloc_onnode");
			exit(1);
		}
		/* copy data */
		memcpy(numa_bcsr->bvalues, bvalues, sizeof(ELEM_TYPE)*nnz);
		bvalues += nnz;
		memcpy(numa_bcsr->bcol_ind, bcol_ind, sizeof(SPM_CRS_IDX_TYPE)*nnz/(r*c));
		//printf("nnz:%lu copy bcol_ind: %lu bcol_ind[0]=%d\n", nnz, sizeof(SPM_CRS_IDX_TYPE)*nnz/(r*c), numa_bcsr->bcol_ind[0]);
		bcol_ind += nnz/(r*c);
		memcpy(numa_bcsr->brow_ptr, brow_ptr, sizeof(SPM_CRS_IDX_TYPE)*(nrows/r+1));
		brow_ptr += nrows/r;
		//printf("thread:%lu cpu:%d bvalues:%p bcol_ind:%p brow_ptr:%p\n", i, spm_thread->cpu, numa_bcsr->bvalues, numa_bcsr->bcol_ind, numa_bcsr->brow_ptr);
		/* fill numa_bcsr */
		numa_bcsr->nr_rows = nrows;
		numa_bcsr->nr_cols = bcsr_mt->bcsr->nr_cols;
		numa_bcsr->nr_nzeros = nnz;
		numa_bcsr->nr_blocks = nnz / (r*c);
		numa_bcsr->nr_brows = nrows / r;
		numa_bcsr->storage = bcsr_mt->bcsr->storage;
		numa_bcsr->br = bcsr_mt->bcsr->br;
		numa_bcsr->bc = bcsr_mt->bcsr->bc;
		/* ugly hack, so that we can use the same spmv routines */
		numa_bcsr->bvalues -= total_nnz;
		numa_bcsr->bcol_ind -= (total_nnz/(r*c));
		numa_bcsr->brow_ptr -= (bcsr_mt->row_start/r);
		/* make the swap */
		bcsr_mt->bcsr = numa_bcsr;
		total_nnz += nnz;
	}

	return spm_mt;
}

static uint64_t bcsr_numa_size(void *spm)
{
	return 0;
}

#define DECLARE_BCSR_NUMA(c, r)  \
static void *\
bcsr_numa_init_mmf_ ## r ## x ## c(char *f,              \
                              unsigned long *rows,       \
                              unsigned long *cols,       \
                              unsigned long *nnz)        \
{                                                        \
   return init_mmf_numa_wrap(f, rows, cols, nnz, r, c);  \
}                                                        \
                                                         \
static void spmv_bcsr_numa_ ## r ## x ## c (void *spm,   \
                                       VECTOR_TYPE *in,  \
                                       VECTOR_TYPE *out) \
{                                                        \
       return bcsr_spmv_wrap(spm, in, out, r, c);        \
}                                                        \
XSPMV_MT_METH_INIT(spmv_bcsr_numa_ ## r ## x ## c,       \
               bcsr_numa_init_mmf_ ## r ## x ## c,       \
               bcsr_numa_size,                           \
               bcsr_destroy,                             \
               sizeof(ELEM_TYPE))                        \

DECLARE_BCSR_NUMA(1,1)
DECLARE_BCSR_NUMA(1,2)
DECLARE_BCSR_NUMA(1,3)
DECLARE_BCSR_NUMA(1,4)
DECLARE_BCSR_NUMA(2,1)
DECLARE_BCSR_NUMA(2,2)
DECLARE_BCSR_NUMA(2,3)
DECLARE_BCSR_NUMA(2,4)
DECLARE_BCSR_NUMA(3,1)
DECLARE_BCSR_NUMA(3,2)
DECLARE_BCSR_NUMA(3,3)
DECLARE_BCSR_NUMA(3,4)
DECLARE_BCSR_NUMA(4,1)
DECLARE_BCSR_NUMA(4,2)
DECLARE_BCSR_NUMA(4,3)
DECLARE_BCSR_NUMA(4,4)

#endif /* SPM_NUMA */

