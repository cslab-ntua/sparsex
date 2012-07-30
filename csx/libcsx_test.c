#include "elmerif.h"
#include "spm_crs.h"
#include "vector.h"
#include <stdio.h>
#include <stdlib.h>

void print_usage(FILE *fs, const char *progname)
{
	fprintf(fs, "%s MM_FILE\n", progname);
}

int main(int argc, char **argv)
{
	if (argc < 2) {
		fprintf(stderr, "%s: too few arguments\n", argv[0]);
		print_usage(stderr, argv[0]);
		exit(1);
	}

	char *mmf_file = argv[1];
	uint64_t nr_rows, nr_cols, nr_nzeros;
	spm_crs32_double_t *csr = (spm_crs32_double_t *)
	    spm_crs32_double_init_mmf(mmf_file,
	                              &nr_rows, &nr_cols, &nr_nzeros, NULL);
	void *tuned = csx_mattune((elmer_index_t *) csr->row_ptr,
	                          (elmer_index_t *) csr->col_ind,
	                          (elmer_value_t *) csr->values,
	                          (elmer_index_t) nr_rows,
	                          (elmer_index_t) nr_cols);

	vector_double_t *x = vector_double_create(nr_cols);
	vector_double_t *y = vector_double_create(nr_rows);
	vector_double_init(x, 1);
	vector_double_init(y, 0);
	csx_matvec(tuned, x->elements, x->size, y->elements, y->size);
	return 0;
}
