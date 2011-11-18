/*
 *  util.h -- Miscellaneous utility stuff (may change in future).
 *
 *  V. K. Karakasis, 2008, bkk@cslab.ece.ntua.gr
 */ 

#ifndef __UTIL_H__
#define __UTIL_H__

#include <limits.h>
#include <stddef.h>

/* Ceiling of two integers */
#define iceil(a, b) ((a) / (b) + ((a) % (b) != 0))

/* Min/Max macros */
#define SPM_MIN(a, b) ((a) < (b) ? (a) : (b))
#define SPM_MAX(a, b) ((a) > (b) ? (a) : (b))

#if SPM_CRS_IDX_TYPE == int
#   define SPM_INDEX_INVALID   INT_MAX
#else
#   define SPM_INDEX_INVALID   LONG_MAX
#endif

/**
 *  Insert a <tt>new</tt> element at position <tt>where</tt> in an array of
 *  <tt>nr_elem</tt> elements, each of size <tt>size</tt>, starting at address
 *  <tt>start</tt>. All elements after <tt>where</tt> are shifted right and
 *  the last element is stored in <tt>overflow</tt>, unless it is <tt>NULL</tt>.
 *
 *  @param start    The starting address of the array.
 *  @param nr_elem  The number of elements in the array.
 *  @param size     The size of each element of the array.
 *  @param where    Where to insert the new element.
 *  @param new      The new element to insert.
 *  @param overflow The rightmost element of the array, that was shifted out
 *                  because of the insertion.
 */
void array_insert(void *start, size_t nr_elem, size_t size,
                  size_t where, const void *new, void *overflow);


#endif  /* __UTIL_H__ */
