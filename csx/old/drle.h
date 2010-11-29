#ifndef __SPM_PATTERNS
#define __SPM_PATTERNS

#include <inttypes.h>

typedef struct {
	uint64_t x,y;
} point_t;

typedef struct {
	point_t  p0;
	uint32_t size;
	uint32_t drle_len;
	enum {DRLE_H=1, DRLE_V, DRLE_D, DRLE_RD} type;
} drle_t;

typedef enum {
	ELEM_SINGLE = 1,
	ELEM_DRLE
} spm_elem_type_t;

typedef struct {
	spm_elem_type_t type;
	union {
		point_t point;   /* ELEM_SIGNLE */
		drle_t *pattern; /* ELEM_DRLE */
	};
} spm_elem_t;

typedef struct {
	uint64_t nrows, ncols, nnz;
	enum {HORIZONTAL=1, VERTICAL, DIAGONAL, RDIAGONAL} type;
	uint64_t elems_nr;
	spm_elem_t *elems;
} spm_t;

#endif /* __SPM_PATTERNS */
