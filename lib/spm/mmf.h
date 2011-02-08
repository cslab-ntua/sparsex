#ifndef __MMF_H__
#define __MMF_H__

#include <stdio.h>
#include <stdint.h>

FILE *mmf_init(char *filename,
               uint64_t *rows, uint64_t *cols,
               uint64_t *nz);

int mmf_get_next(FILE *mmf,
                 uint64_t *row, uint64_t *col,
                 double *val);

int mmf_get_next_vstr(FILE *mmf,
                      uint64_t *row, uint64_t *col,
                      char **val);
#endif
