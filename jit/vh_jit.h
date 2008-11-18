#ifndef VH_JIT__
#define VH_JIT__

#include "huffman.h"

void *vhjit_init(huff_node_t *tree, unsigned long bits, hj_out_t type);
void *vhjit_get_spmvfn(void *hj_wrap, int ci_bits);
void vhjit_shut(void *hj_wrap);

#endif /* VH_JIT__ */
