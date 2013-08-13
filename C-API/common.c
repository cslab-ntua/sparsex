#include "common.h"

void *malloc_internal(size_t x, const char *sourcefile, unsigned long lineno,
                      const char *function)
{
    void *ret;
    ret = malloc(x);
    if (!ret) {
        err_handle(LIBCSX_ERR_MEM_ALLOC, sourcefile, lineno, function, NULL);
        exit(1);
    }
    return ret;
}

void free_internal(void *ptr, const char *sourcefile, unsigned long lineno,
                   const char *function)
{
    if (!ptr) {
        err_handle(LIBCSX_ERR_MEM_FREE, sourcefile, lineno, function, NULL);
        exit(1);
    }
    free(ptr);
    ptr = NULL;
}
