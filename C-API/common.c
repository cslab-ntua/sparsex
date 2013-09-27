#include "common.h"

void libcsx_init()
{
#ifdef LIBCSX_LOG_FILE
    err_set_logfile("libcsx_logfile");
#else
    err_set_logfile(NULL);    
#endif
}

void libcsx_close()
{
#ifdef LIBCSX_LOG_FILE
    err_close_logfile();
#endif
}

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
