/*
 *  util.c -- Implementation of the spm/util.h interface.
 *
 *  V. K. Karakasis, 2008, bkk@cslab.ece.ntua.gr
 */ 

#include <string.h>
#include "util.h"

void
array_insert(void *start, size_t nr_elem, size_t size,
             size_t where, const void *new, void *overflow)
{
    size_t  i;

    if (overflow)
        memcpy(overflow, start + size*(nr_elem - 1), size);

    for (i = nr_elem - 1; i > where; i--)
        memcpy(start + size * i, start + size*(i - 1), size);

    memcpy(start + size*where, new, size);
    return;
}
