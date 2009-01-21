#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <string.h>

#include "spm_crs.h"


int main(int argc, char **argv)
{
	if (argc < 2){
		fprintf(stderr, "Usage: %s <mmf_file>\n", argv[0]);
		exit(1);
	}

	char *mmf_file = argv[1];
	uint64_t nrows0, ncols0, nnz0;
	uint64_t nrows1, ncols1, nnz1;

	spm_crs32_double_t *crs0, *crs1;
	crs0 = spm_crs32_double_init_mmf(mmf_file, &nrows0, &ncols0, &nnz0);
	crs1 = spm_crs32_double_init_mmf(mmf_file, &nrows1, &ncols1, &nnz1);

	assert(nrows0=nrows1);
	assert(ncols0=ncols1);
	assert(nnz0=nnz1);

	assert(memcmp(crs0->row_ptr, crs1->row_ptr, nrows0+1) == 0);
	assert(memcmp(crs0->col_ind, crs1->col_ind, nnz0) == 0);
	assert(memcmp(crs0->values, crs1->values, nnz0) == 0);

	return 0;
}
