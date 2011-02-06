#ifndef __EXT_PROG_H__
#define __EXT_PROG_H__

#include <stdio.h>

struct extp_s;
typedef struct extp_s extp_t;

extp_t *extp_open(char *ext_prog, FILE **in, FILE **out);
int extp_close(extp_t *extp, FILE *in, FILE *out);

#endif /* __EXT_PROG_H__ */
