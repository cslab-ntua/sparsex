#ifndef __MMF_H__
#define __MMF_H__

#include <stdio.h>

FILE *mmf_init(char *filename, 
               unsigned long *rows, unsigned long *cols, 
               unsigned long *nz);

int mmf_get_next(FILE *mmf, 
                 unsigned long *row, unsigned long *col, 
                 double *val);

int mmf_get_next_vstr(FILE *mmf, 
                      unsigned long *row, unsigned long *col, 
                      char **val);
#endif
