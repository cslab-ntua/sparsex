/*
 *  spm/bitstr.h -- A simple and portable bit string implementation.
 *
 *  For a more sophisticated and efficient implementation, you should look at
 *  specific processor instructions that test/set a bit in memory.
 *
 *  V. K. Karakasis, 2008, bkk@cslab.ece.ntua.gr
 */ 
#ifndef __BITSTR_H__
#define __BITSTR_H__

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

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

/*
 *  FIXME: *alloc can fail silently.
 */ 
inline static bitstr_t *
bs_create(unsigned long nr_bits)
{
    bitstr_t        *new = malloc(sizeof(*new));
    unsigned long   nr_words = iceil(nr_bits, LONG_BIT);

    new->base_addr = calloc(nr_words, sizeof(*new->base_addr));
    new->nr_bits   = nr_words * LONG_BIT;
    return new;
}

inline static void
bs_delete(bitstr_t *bs)
{
    if (!bs)
        return;

    free(bs->base_addr);
    free(bs);
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

#if LONG_BIT_DEFINED
#   undef  LONG_BIT     /* We defined LONG_BIT; undefine it. */
#endif

#undef  BS_WORD     /* That was in internal macro; not needed any more */
#endif              /* __BITSTR_H__ */
