#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "beanbag.h"
#include "method.h"
#include "spm_crsr.h"
#include "spm_parse.h"
#include "mmf.h"
#include "tsc.h"
//#include "prefetch.h"

spm_crsr_t *spm_crsr_create(unsigned long nrows, 
                            unsigned long ncols,
                            unsigned long nz ) 
{
	spm_crsr_t *matrix;
	matrix = malloc(sizeof(spm_crsr_t));
	if ( !matrix ){
		perror("malloc");
		exit(1);
	}
	//printf("values: %lu\n", (unsigned long)&matrix->values);
	matrix->values = malloc(sizeof(elem_t)*nz);
	if ( !matrix->values ){
		perror("malloc");
		exit(1);
	}
	/*
	*/
	matrix->rows = calloc(nrows, sizeof(spm_crsr_row_t));
	if ( !matrix->rows ) {
		perror("malloc");
		exit(1);
	}
		 
	matrix->nrows = nrows;
	matrix->ncols = ncols;
	matrix->nnz = nz;
	/*
	*/
	return matrix;
}

spm_crsr2_t *spm_crsr2_create(unsigned long nrows,
                              unsigned long ncols,
                              unsigned long nz)
{
	spm_crsr2_t *matrix;
	matrix = malloc(sizeof(spm_crsr2_t));
	if ( !matrix ){
		perror("malloc");
		exit(1);
	}
	//printf("values: %lu\n", (unsigned long)&matrix->values);
	matrix->values = malloc(sizeof(elem_t)*nz);
	if ( !matrix->values ){
		perror("malloc");
		exit(1);
	}

	matrix->nrows = nrows;
	matrix->ncols = ncols;
	matrix->nnz = nz;

	return matrix;
}
                              

void spm_crsr_destroy(spm_crsr_t *matrix) 
{
     free(matrix->values);
     free(matrix->rows);
     free(matrix);
}

void spm_crsr_print(spm_crsr_t *crsr)
{
	spm_crsr_row_t *row;
	spm_crsr_unit_t *unit;

	unsigned long i, j, k;
	for ( i=0; i < crsr->nrows; i++ ){
		row = &crsr->rows[i];
		printf("row: %lu (%p)\n", i, row);
		for (j=0; j < row->units_nr; j++){
			unit = &row->units[j];
			printf(" unit: %lu (%p)\n", j, unit);
			for ( k=0; k < unit->len; k++){
				printf("elem: %lu %lf colind:%lu\n", k, crsr->values[k+unit->val_indx], k+unit->col_indx);
			}
		}
	}
	printf("\n");
}

void spm_crsr2_print(spm_crsr2_t *crsr2)
{
	spm_crsr2_row_t *row;
	spm_crsr_unit_t *unit;

	unsigned long i, j, k;
	for ( i=0; i < crsr2->nrows; i++ ){
		row = &crsr2->rows[i];
		printf("row: %lu (%p) %lu \n", i, row, row->row_indx);
		for (j=0; j < row->units_nr; j++){
			unit = &row->units[j];
			printf("\tunit: %lu (%p) len:%lu col_ind:%lu val_ind:%lu\n", 
			       j, unit, unit->len, unit->col_indx, unit->val_indx);
			#if 0
			for ( k=0; k < unit->len; k++){
				printf("elem: %lu %lf colind:%lu\n", k, crsr2->values[k+unit->val_indx], k+unit->col_indx);
			}
			#endif
		}
	}
	printf("\n");
}

static inline void *alloc_elems_aligned(unsigned long elems)
{
	void *data;
	int err;
	
	err = posix_memalign(&data, 16, elems*sizeof(elem_t));
	if (err){
		fprintf(stderr, "posix_memalign error: %d", err);
		exit(1);
	}

	return data;
}

#define BEANBAG_ALLOC_UNITS 4096
#define BEANBAG_ALLOC_ELEMS 4096
/*
 * CRSR Format Functions
 */
struct crsr_state_s {
	spm_crsr_t *crsr;
	beanbag_t *units_bb;
	unsigned long val_indx;
};
typedef struct crsr_state_s crsr_state_t;


static void spm_crsr_init(parse_state_t *parse_state, void **crsr_state_ptr,
                          unsigned long rows_nr, unsigned long cols_nr,
                          unsigned long nz_nr)
{
	spm_crsr_t *crsr;
	crsr_state_t *crsr_state;

	crsr = spm_crsr_create(rows_nr, cols_nr, nz_nr);
	crsr_state = malloc(sizeof(crsr_state_t));
	if ( !crsr || !crsr_state ){
		perror("spm_crsr_init: malloc");
		exit(1);
	}
	
	crsr_state->crsr = crsr;
	crsr_state->val_indx = 0;
	crsr_state->units_bb = beanbag_create(sizeof(spm_crsr_unit_t), BEANBAG_ALLOC_UNITS);

	*crsr_state_ptr  = crsr_state;
}

static void *spm_crsr_finalize(parse_state_t *parse_state, void *crsr_st)
{
	crsr_state_t *crsr_state = crsr_st;
	spm_crsr_t *crsr = crsr_state->crsr;

	crsr->nnz = crsr_state->val_indx;
	//crsr->values = realloc(crsr->values, crsr->nnz);

	free(crsr_state);
	return crsr;
}

static void spm_crsr_finalize_unit(parse_state_t *parse_state, void *crsr_st)
{
	unsigned long l;
	crsr_state_t *crsr_state = crsr_st;
	spm_crsr_t *crsr = crsr_state->crsr;
	spm_crsr_unit_t *u = beanbag_alloc(crsr_state->units_bb);

	u->col_indx = parse_state->cur_unit_col;
	u->len = parse_state->cur_values_idx;
	u->val_indx = crsr_state->val_indx;

	for (l=0; l < parse_state->cur_values_idx; l++){
		crsr->values[u->val_indx + l] = parse_state->cur_values[l];
	}
	crsr_state->val_indx += parse_state->cur_values_idx;
}

static void spm_crsr_finalize_row(parse_state_t *parse_state, void *crsr_st)
{
	crsr_state_t *crsr_state = crsr_st;
	spm_crsr_row_t *r;

	r = &(crsr_state->crsr->rows[parse_state->cur_row]);
	r->units_nr = beanbag_size(crsr_state->units_bb);
	r->units = beanbag_destroy(crsr_state->units_bb);
	crsr_state->units_bb = beanbag_create(sizeof(spm_crsr_unit_t), BEANBAG_ALLOC_UNITS);
}

static spm_parser_t crsr_parser = {
	.initialize = spm_crsr_init,
	.finalize = spm_crsr_finalize,
	.finalize_unit = spm_crsr_finalize_unit,
	.finalize_row = spm_crsr_finalize_row
};

/*
 * CRSR2 Format Functions
 */
struct crsr2_state_s {
	spm_crsr2_t  *crsr2;
	beanbag_t *units_bb, *rows_bb;
	unsigned long val_indx;
	unsigned long row_indx;
};
typedef struct crsr2_state_s crsr2_state_t;

static void spm_crsr2_init(parse_state_t *parse_state, void **crsr2_state_ptr,
                           unsigned long rows_nr, unsigned long cols_nr, 
                           unsigned long nz_nr)
{
	spm_crsr2_t *crsr2;
	crsr2_state_t *crsr2_state;

	crsr2 = spm_crsr2_create(rows_nr, cols_nr, nz_nr);
	crsr2_state = malloc(sizeof(crsr2_state_t));
	if ( !crsr2 || !crsr2_state ){
		perror("spm_crsr2_init: malloc");
		exit(1);
	}

	crsr2_state->crsr2 = crsr2;
	crsr2_state->units_bb = beanbag_create(sizeof(spm_crsr_unit_t), BEANBAG_ALLOC_UNITS);
	crsr2_state->rows_bb = beanbag_create(sizeof(spm_crsr2_row_t), crsr2->nrows);
	crsr2_state->val_indx = 0;
	crsr2_state->row_indx = 0;
	*crsr2_state_ptr = crsr2_state;
}

static void *spm_crsr2_finalize(parse_state_t *parse_state, void *crsr2_st)
{
	crsr2_state_t *crsr2_state = crsr2_st;
	spm_crsr2_t *crsr2 = crsr2_state->crsr2;

	crsr2->nnz = crsr2_state->val_indx;
	//crsr2->values = realloc(crsr2->values, crsr2->nnz);

	crsr2->nrows = beanbag_size(crsr2_state->rows_bb);
	crsr2->rows = beanbag_destroy(crsr2_state->rows_bb);

	free(crsr2_state);
	return crsr2;
}

static void spm_crsr2_finalize_unit(parse_state_t *parse_state, void *crsr2_st)
{
	unsigned long l;
	crsr2_state_t *crsr2_state = crsr2_st;
	spm_crsr2_t *crsr2 = crsr2_state->crsr2;
	spm_crsr_unit_t *u = beanbag_alloc(crsr2_state->units_bb);

	u->col_indx = parse_state->cur_unit_col;
	u->len = parse_state->cur_values_idx;
	u->val_indx = crsr2_state->val_indx;

	for (l=0; l < parse_state->cur_values_idx; l++){
		crsr2->values[u->val_indx + l] = parse_state->cur_values[l];
	}
	crsr2_state->val_indx += parse_state->cur_values_idx;
	parse_state->cur_values_idx = 0;
}

static void spm_crsr2_finalize_row(parse_state_t *parse_state, void *crsr2_st)
{
	crsr2_state_t *crsr2_state = crsr2_st;
	unsigned long units_nr;
	units_nr = beanbag_size(crsr2_state->units_bb);
	if ( units_nr ){
		spm_crsr2_row_t *r = beanbag_alloc(crsr2_state->rows_bb);
		r->row_indx = crsr2_state->row_indx;
		r->units_nr = units_nr;
		r->units = beanbag_destroy(crsr2_state->units_bb);
		crsr2_state->units_bb = beanbag_create(sizeof(spm_crsr_unit_t), BEANBAG_ALLOC_UNITS);
	}
	crsr2_state->row_indx++;
}

static spm_parser_t crsr2_parser = {
	.initialize = spm_crsr2_init,
	.finalize = spm_crsr2_finalize,
	.finalize_unit = spm_crsr2_finalize_unit,
	.finalize_row = spm_crsr2_finalize_row
};


/* XXX: hybrid is broken ATM */
#if 0


/*
 * Hybrid Format Functions
 */

struct hybrid_state_s {
	crsr_state_t *crsr;
	crs_state_t  *crs;
	spm_hybrid_t *hybrid;
};
typedef struct hybrid_state_s hybrid_state_t;

static void spm_hybrid_initialize(parse_state_t *parse_state, void **hybrid_state_ptr,
                                  unsigned long rows_nr, unsigned long cols_nr,
                                  unsigned long nz_nr)
{
	hybrid_state_t *hybrid_state;
	spm_hybrid_t *hybrid;
	void *crs, *crsr;

	hybrid_state = malloc(sizeof(hybrid_state_t));
	hybrid = malloc(sizeof(spm_hybrid_t));
	if (!hybrid_state || !hybrid){
		perror("spm_hybrid_initialize: malloc");
		exit(1);
	}

	hybrid_state->hybrid = hybrid;
	spm_crs_initialize(parse_state, &crs, rows_nr, cols_nr, nz_nr);
	hybrid_state->crs = crs;
	spm_crsr_init(parse_state, &crsr, rows_nr, cols_nr, nz_nr);
	hybrid_state->crsr = crsr;

	*hybrid_state_ptr = hybrid_state;
}

static void *spm_hybrid_finalize(parse_state_t *parse_state, void *hybrid_st)
{
	hybrid_state_t *hybrid_state = hybrid_st;
	spm_hybrid_t *hybrid = hybrid_state->hybrid;

	hybrid->crsr = spm_crsr_finalize(parse_state, hybrid_state->crsr);
	hybrid->crs = spm_crs_finalize(parse_state, hybrid_state->crs);

	free(hybrid_state);
	return hybrid;
}

static void spm_hybrid_finalize_unit(parse_state_t *parse_state, void *hybrid_st)
{
	hybrid_state_t *hybrid_state = hybrid_st;
	if ( parse_state->cur_values_idx < 4){
		spm_crs_finalize_unit(parse_state, hybrid_state->crs);
	} else {
		spm_crsr_finalize_unit(parse_state, hybrid_state->crsr);
	}
}

static void spm_hybrid_finalize_row(parse_state_t *parse_state, void *hybrid_st)
{
	hybrid_state_t *hybrid_state = hybrid_st;
	spm_crs_finalize_row(parse_state, hybrid_state->crs);
	spm_crsr_finalize_row(parse_state, hybrid_state->crsr);
}


static spm_parser_t hybrid_parser = {
	.initialize = spm_hybrid_initialize,
	.finalize = spm_hybrid_finalize,
	.finalize_unit = spm_hybrid_finalize_unit,
	.finalize_row = spm_hybrid_finalize_row
};

/*
 * Hybrid2 Format Functions
 */

struct hybrid2_state_s {
	crsr_state_t *crsr2;
	crs_state_t  *crs;
	spm_hybrid2_t *hybrid2;
};
typedef struct hybrid2_state_s hybrid2_state_t;

static void spm_hybrid2_initialize(parse_state_t *parse_state, void **hybrid2_state_ptr,
                                  unsigned long rows_nr, unsigned long cols_nr,
                                  unsigned long nz_nr)
{
	hybrid2_state_t *hybrid2_state;
	spm_hybrid2_t *hybrid2;
	void *crs, *crsr2;

	hybrid2_state = malloc(sizeof(hybrid2_state_t));
	hybrid2 = malloc(sizeof(spm_hybrid2_t));
	if (!hybrid2_state || !hybrid2){
		perror("spm_hybrid2_initialize: malloc");
		exit(1);
	}

	hybrid2_state->hybrid2 = hybrid2;
	spm_crs_initialize(parse_state, &crs, rows_nr, cols_nr, nz_nr);
	hybrid2_state->crs = crs;
	spm_crsr2_init(parse_state, &crsr2, rows_nr, cols_nr, nz_nr);
	hybrid2_state->crsr2 = crsr2;

	*hybrid2_state_ptr = hybrid2_state;
}

static void *spm_hybrid2_finalize(parse_state_t *parse_state, void *hybrid2_st)
{
	hybrid2_state_t *hybrid2_state = hybrid2_st;
	spm_hybrid2_t *hybrid2 = hybrid2_state->hybrid2;

	hybrid2->crsr2 = spm_crsr2_finalize(parse_state, hybrid2_state->crsr2);
	hybrid2->crs = spm_crs_finalize(parse_state, hybrid2_state->crs);

	free(hybrid2_state);
	return hybrid2;
}

static void spm_hybrid2_finalize_unit(parse_state_t *parse_state, void *hybrid2_st)
{
	hybrid2_state_t *hybrid2_state = hybrid2_st;
	if ( parse_state->cur_values_idx < 4){
		spm_crs_finalize_unit(parse_state, hybrid2_state->crs);
	} else {
		spm_crsr2_finalize_unit(parse_state, hybrid2_state->crsr2);
	}
}

static void spm_hybrid2_finalize_row(parse_state_t *parse_state, void *hybrid2_st)
{
	hybrid2_state_t *hybrid2_state = hybrid2_st;
	spm_crs_finalize_row(parse_state, hybrid2_state->crs);
	spm_crsr2_finalize_row(parse_state, hybrid2_state->crsr2);
}

static spm_parser_t hybrid2_parser = {
	.initialize = spm_hybrid2_initialize,
	.finalize = spm_hybrid2_finalize,
	.finalize_unit = spm_hybrid2_finalize_unit,
	.finalize_row = spm_hybrid2_finalize_row
};


void print_spm_hybrid_stats(spm_hybrid_t *hybrid, 
                            unsigned long rows_nr, unsigned long cols_nr,
			    unsigned long nz_nr)
{
	printf("mat : nz:%lu rows:%lu cols:%lu\n", nz_nr, rows_nr, cols_nr);
	printf("crsr: nz:%lu \n", hybrid->crsr->nnz);
}


spm_hybrid_t *spm_hybrid_init_mmf(char *mmf_file,
                                  unsigned long *rows_nr, unsigned long *cols_nr,
                                  unsigned long *nz_nr)
{
	int ret;
	spm_hybrid_t *hybrid;
	FILE *f;

	f = mmf_init(mmf_file, rows_nr, cols_nr, nz_nr);
	hybrid = spm_parse_seq(f, *rows_nr, *cols_nr, *nz_nr, &hybrid_parser);
	fclose(f);
	print_spm_hybrid_stats(hybrid, *rows_nr, *cols_nr, *nz_nr);

	return hybrid;
}

void spm_hybrid_stats(char *mmf_file)
{
	unsigned long rows_nr, cols_nr, nz_nr;
	spm_hybrid_t *hybrid = spm_hybrid_init_mmf(mmf_file, &rows_nr, &cols_nr, &nz_nr);
	print_spm_hybrid_stats(hybrid, rows_nr, cols_nr, nz_nr);
}

void print_spm_hybrid2_stats(spm_hybrid2_t *hybrid2, 
                            unsigned long rows_nr, unsigned long cols_nr,
			    unsigned long nz_nr)
{
	printf("mat : nz:%lu rows:%lu cols:%lu\n", nz_nr, rows_nr, cols_nr);
	printf("crsr2: nz:%lu \n", hybrid2->crsr2->nnz);
	spm_crsr2_print(hybrid2->crsr2);
}

spm_hybrid2_t *spm_hybrid2_init_mmf(char *mmf_file,
                                  unsigned long *rows_nr, unsigned long *cols_nr,
                                  unsigned long *nz_nr)
{
	int ret;
	spm_hybrid2_t *hybrid2;
	FILE *f;

	f = mmf_init(mmf_file, rows_nr, cols_nr, nz_nr);
	hybrid2 = spm_parse_seq(f, *rows_nr, *cols_nr, *nz_nr, &hybrid2_parser);
	fclose(f);
	print_spm_hybrid2_stats(hybrid2, *rows_nr, *cols_nr, *nz_nr);

	return hybrid2;
}



void spm_crsr2_multiply(spm_crsr2_t *crsr2, vector_t *in, vector_t *out)
{
	elem_t *x = in->elements, *y = out->elements, *values = crsr2->values;
	const unsigned long nrows = crsr2->nrows;
	unsigned long r, u, i;
	spm_crsr2_row_t *row = crsr2->rows;

	for (r=0; r< nrows; r++){
		const unsigned long units_nr = row->units_nr;
		spm_crsr_unit_t *unit = row->units;
		register elem_t yr=0;
		for (u=0; u< units_nr; u++){
			const unsigned long unit_len = unit->len;
			const elem_t *myx = &x[unit->col_indx];
			const elem_t *myv = &values[unit->val_indx];
			for ( i=0; i<unit_len; i++){
				yr += myx[i]*myv[i];
			}
			unit++;
		}
		y[row->row_indx] += yr;
		row++;
	}
}

#if 0
static tsc_t tsc_crs, tsc_crsr;
void __attribute__((constructor)) __tinit(void)
{
	tsc_init(&tsc_crs);
	tsc_init(&tsc_crsr);
}
void __attribute__((destructor)) __tdest(void)
{
	printf("crs:\t");
	tsc_report(&tsc_crs);
	printf("crsr(2):\t");
	tsc_report(&tsc_crsr);
}
#endif

void spm_crs_multiply_safe(spm_t *matrix, vector_t *in, vector_t *out);
void spm_hybrid_multiply(spm_hybrid_t *hybrid, vector_t *in, vector_t *out)
{
	//tsc_start(&tsc_crs);
	spm_crs_multiply_safe((spm_t *)hybrid->crs, in, out);
	//tsc_pause(&tsc_crs);

	//tsc_start(&tsc_crsr);
	spm_crsr_multiply(hybrid->crsr, in, out);
	//tsc_pause(&tsc_crsr);

}
METHOD_INIT(spm_hybrid_multiply, spm_hybrid_init_mmf)

void spm_hybrid2_multiply(spm_hybrid2_t *hybrid2, vector_t *in, vector_t *out)
{
	//tsc_start(&tsc_crsr);
	spm_crsr2_multiply(hybrid2->crsr2, in, out);
	//tsc_pause(&tsc_crsr);

	//tsc_start(&tsc_crs);
	spm_crs_multiply_safe((spm_t *)hybrid2->crs, in, out);
	//tsc_pause(&tsc_crs);
}
METHOD_INIT(spm_hybrid2_multiply, spm_hybrid2_init_mmf)

void spm_hybrid_multiply2(spm_hybrid_t *hybrid, vector_t *in, vector_t *out)
{
	elem_t *x = in->elements, *y = out->elements;
	spm_crs_t *crs = hybrid->crs;
	spm_crsr_t *crsr = hybrid->crsr;
	elem_t *crs_values = hybrid->crs->values;
	elem_t *crsr_values = hybrid->crsr->values;
	index_t *row_ptr = hybrid->crs->row_ptr, *col_ind = hybrid->crs->col_ind;
	spm_crsr_row_t *row = crsr->rows;
	unsigned long nrows = crsr->nrows;
	unsigned long r, u, i;
	register elem_t yr;

	for ( r=0; r< nrows; r++){
		const unsigned long units_nr = row->units_nr;
		spm_crsr_unit_t *unit = row->units;
		yr=0;
		for ( u=0 ; u < units_nr; u++){
			const unsigned long unit_len = unit->len;
			const elem_t *myx = &x[unit->col_indx];
			const elem_t *myv = &crsr_values[unit->val_indx];
			for ( i=0; i<unit_len; i++){
				yr += myx[i]*myv[i];
			}
			unit++;
		}
		for (u=row_ptr[r]; u<row_ptr[r+1]; u++){
			yr += (crs_values[u] * x[col_ind[u]]);
		}
		row++;
		y[r] = yr;
	}
}
METHOD_INIT(spm_hybrid_multiply2, spm_hybrid_init_mmf)
#endif /* broken hybrid */


spm_crsr_t *spm_crsr_init_mmf(char *mmf_file,
                              unsigned long *rows_nr, unsigned long *cols_nr,
                              unsigned long *nz_nr)
{
	int ret;
	spm_crsr_t *crsr;
	FILE *f;
	
	f = mmf_init(mmf_file, rows_nr, cols_nr, nz_nr);
	crsr = spm_parse_seq(f, *rows_nr, *cols_nr, *nz_nr, &crsr_parser);
	fclose(f);

	return crsr;
}


void spm_crsr_multiply(spm_crsr_t *crsr, vector_t *in, vector_t *out)
{
	elem_t *x = in->elements, *y = out->elements, *values = crsr->values;
	const unsigned long nrows = crsr->nrows;
	unsigned long r, u, i;
	spm_crsr_row_t *row = crsr->rows;
	register elem_t yr;

	for (r=0; r< nrows; r++){
		const unsigned long units_nr = row->units_nr;
		spm_crsr_unit_t *unit = row->units;
		yr=0;
		//printf("-----> row=%lu\n", r);
		for (u=0; u< units_nr; u++){
			const unsigned long unit_len = unit->len;
			const elem_t *myx = &x[unit->col_indx];
			const elem_t *myv = &values[unit->val_indx];
			for ( i=0; i<unit_len; i++){
				//printf("\t%f += (%f*%f) = ", yr, myx[i], myv[i]);
				yr += myx[i]*myv[i];
				//printf("%f\n", yr);
			}
			unit++;
		}
		row++;
		y[r] += yr;
	}
}
METHOD_INIT(spm_crsr_multiply, spm_crsr_init_mmf)

#ifdef SPM_CRSR_TEST

static void spm_crsr_test(FILE *f, spm_crsr_t *crsr){
	unsigned long r, u, i;
	unsigned long ro, co;
	elem_t eo;
	for ( r=0; r < crsr->nrows; r++){
		spm_crsr_row_t *row = &crsr->rows[r];
		for ( u=0; u < row->units_nr; u++ ){
			spm_crsr_unit_t *unit = &row->units[u];
			for (i=0 ; i < unit->len ; i++){
				elem_t e = crsr->values[unit->val_indx + i];
				mmf_get_next(f, &ro, &co, &eo);
				if (!((r==ro) && ((unit->col_indx + i) == co) && (e == eo))){
					fprintf(stderr, "\nvalidation error:\n");	
					fprintf(stderr, "  r=%lu ro=%lu\n", r, ro);	
					fprintf(stderr, "  c=%lu co=%lu\n", unit->col_indx + i, co);	
					fprintf(stderr, "  e=%lf eo=%lf\n", e, eo);	
					exit(1);
				}
			}
		}
	}
}

int main(int argc, char **argv)
{
	spm_crsr_t *crsr;
	unsigned long rows, cols, nz;
	unsigned long nrows, ncols, nnz;
	FILE *f;
	crsr = spm_crsr_init_mmf(argv[1], &nrows, &ncols, &nnz);
	///printf("veryfing ...");
	f = mmf_init(argv[1], &rows, &cols, &nz);
	spm_crsr_test(f, crsr);
	//printf("... verified\n");
	
	return 0;
}
#endif

