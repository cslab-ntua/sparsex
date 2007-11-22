#include <stdlib.h>
#include <stdio.h>

#ifdef NUMA
#include <numa.h>
#endif

#include "matrix.h"

/*
 * FIXME: a number of functions have not been converted to the new
 *        dynamic element type manipulation
 */

matrix_t *matrix_create(unsigned long nrows, unsigned long ncols, unsigned elem_type) 
{
	unsigned long i;
	matrix_t *matrix;

	matrix = malloc(sizeof(matrix_t));
	if ( !matrix ){
		perror("malloc: matrix");
		exit(1);
	}
	matrix->nrows = nrows;
	matrix->ncols = ncols;
	matrix->elem_type = elem_type;

	matrix->elements = malloc(elem_size(elem_type)*nrows*ncols);
	matrix->rows = malloc(sizeof(void *)*nrows);
	if ( !matrix->elements || matrix->rows ){
		perror("malloc: matrix->elements || matrix->rows");
		exit(1);
	}

	for (i=0; i<nrows; i++) {
		matrix->rows[i] = matrix->elements + i*ncols*elem_size(elem_type);
	}

	return matrix;
}

void matrix_destroy(matrix_t *matrix) 
{
	free(matrix->elements);
	free(matrix->rows);
	free(matrix);
}


vector_t *vector_create(unsigned long size, unsigned elem_type) 
{
	vector_t *vector = malloc(sizeof(vector_t));
	vector->size = size;
	vector->elements = malloc(elem_size(elem_type)*size);
	return vector;
}

void vector_destroy(vector_t *vector) 
{
	free(vector->elements);
	free(vector);
}

#ifdef NUMA
vector_t *numa_vector_create(unsigned long size, int memnode) 
{
	vector_t *vector; 
	vector = numa_alloc_onnode_safe(sizeof(vector_t), memnode_number);
	vector->size = size;
	vector->elements = numa_alloc_onnode_safe(elem_size(elem_type)*size, memnode);
	return vector;
}

void numa_vector_destroy(vector_t *vector) 
{
	numa_free(vector->elements, vector->size*elem_size(elem_type));
	numa_free(vector, sizeof(vector_t));
}
#endif


void matrix_init_float(matrix_t *matrix, float val) 
{
	 unsigned long r, c;
	 float **rows = (float **)matrix->rows;
	 for (r = 0; r < matrix->nrows; r++) {
		  for (c = 0; c < matrix->ncols; c++) {
			   rows[r][c] = val;
		  }
	 }
	 return;
}

void matrix_init_double(matrix_t *matrix, double val) 
{
	 unsigned long r, c;
	 double **rows = (double **)matrix->rows;
	 for (r = 0; r < matrix->nrows; r++) {
		  for (c = 0; c < matrix->ncols; c++) {
			   rows[r][c] = val;
		  }
	 }
	 return;
}

void matrix_init(matrix_t *matrix, double val)
{
	switch (matrix->elem_type){
		case ELEM_FLOAT:
		matrix_init_float(matrix, (float)val);
		break;

		case ELEM_DOUBLE:
		matrix_init_double(matrix, val);
		break;

		default:
		fprintf(stderr, "%s: uknown type: %u\n", __FUNCTION__, matrix->elem_type);
		exit(1);
	}
}

void vector_init_double(vector_t *vector, double val) 
{
	 unsigned long i;

	 double *elements = vector->elements;
	 for(i=0; i<vector->size; i++)
		elements[i] = val;

	 return;
}

void vector_init_float(vector_t *vector, float val) 
{
	 unsigned long i;

	 double *elements = vector->elements;
	 for(i=0; i<vector->size; i++)
		elements[i] = val;

	 return;
}

void vector_init(vector_t *vector, double val)
{
	switch (vector->elem_type){
		case ELEM_FLOAT:
		vector_init_float(vector, (float)val);
		break;

		case ELEM_DOUBLE:
		vector_init_double(vector, val);
		break;

		default:
		fprintf(stderr, "%s: uknown type: %u\n", __FUNCTION__, vector->elem_type);
		exit(1);
	}
}

void vector_init_rand_range(vector_t *vector, double max, double min)
{
	unsigned long i;
	switch (vector->elem_type) {

		case ELEM_FLOAT:
		{
			float *elements = vector->elements;
			float val;
			float fmax = (float)max, fmin = (float)min;
			for ( i=0; i < vector->size; i++){
				val = ((float)(rand()+i)/((float)RAND_MAX+1));
				elements[i] = fmin + val*(fmax-fmin);
			}
		}
		break;

		case ELEM_DOUBLE:
		{
			double *elements = vector->elements;
			double val;
			for ( i=0; i < vector->size; i++){
				val = ((double)(rand()+i)/((double)RAND_MAX+1));
				elements[i] = min + val*(max-min);
			}
		}
		break;

		default:
		fprintf(stderr, "%s: uknown type: %u\n", __FUNCTION__, vector->elem_type);
		exit(1);
	}
}

#if 0
/*
 * Create unitary matrix
 * (matrix is supposed to be square)
 */ 
void matrix_init_unit(matrix_t *matrix) 
{
	 unsigned long r;

	 matrix_init(matrix, 0);

	 for (r = 0; r < matrix->nrows; r++) {
		  matrix->rows[r][r] = 1;
	 }

	 return;
}
#endif


#if 0
/*
 * Initialize subsequent elements with increasing values
 * modulo 'max_val', starting from 0
 *
 */ 
void matrix_init_inc(matrix_t *matrix, int max_val) 
{
	 unsigned long r, c;
	 elem_t val = 0;
	 for (r = 0; r < matrix->nrows; r++) {
		  for (c = 0; c < matrix->ncols; c++) {
			   matrix->rows[r][c] = (int)val % max_val;
			   val++;
		  }
	 }
	 return;
}
#endif


#if 0
/*
 * Initialize matrix elements pseudo-randomly,
 * modulo 'max_val'
 *
 */ 
void matrix_init_rand(matrix_t *matrix, int max_val) 
{
	 unsigned long r, c;
	 elem_t val;

	 for (r = 0; r < matrix->nrows; r++) {
		  for (c = 0; c < matrix->ncols; c++) {
			   val = rand() % max_val;
			   matrix->rows[r][c] = val;
		  }
	 }
	 return;
}

void vector_init_rand(vector_t *vector, int max_val) 
{
	 unsigned long i;
	 elem_t val;

	 for(i=0; i< vector->size; i++) {
		 val = rand() % max_val;
		 vector->elements[i] = val;
	 }
	 
	 return;
}
#endif




#if 0
void matrix_init_rand_range(matrix_t *matrix, elem_t max, elem_t min) 
{
	 unsigned long r, c;
	 elem_t val;

	 for (r = 0; r < matrix->nrows; r++) {
		  for (c = 0; c < matrix->ncols; c++) {
		val = ((elem_t)rand()/((elem_t)RAND_MAX+1));
			   matrix->rows[r][c] = val;
		  }
	 }
	 return;
}
#endif

#if 0
/*
 * LAYOUT ROUTINES
 *
 */ 

/*
 *  Convert matrix m_real, which is stored row-wise, into 
 *  block ZZ layout
 */
void matrix_row2block_ZZ(matrix_t *m_real, matrix_t *m, unsigned int nrows_tile, unsigned int ncols_tile)
{
	int rr, cc, r, c, count;
	unsigned long nrows = m_real->nrows;
	unsigned long ncols = m_real->ncols;

	count = 0;
	for (rr=0; rr < nrows; rr += nrows_tile)
		for (cc=0; cc < ncols; cc += ncols_tile)
			for (r = rr; r < rr + nrows_tile; r++)
				for (c = cc; c < cc + ncols_tile ; c++ )
					m->elements[count++] = m_real->rows[r][c];
}


/*
 *   Convert a table which is stored in block ZZ layout to row-major storage
 */
void matrix_block2row_ZZ(matrix_t *m, matrix_t *m_real, unsigned int nrows_tile, unsigned int ncols_tile)
{
	int rr, cc, r, c, count;
	unsigned long nrows = m_real->nrows;
	unsigned long ncols = m_real->ncols;

	count = 0;
	for (rr=0; rr < nrows; rr += nrows_tile)
		for (cc=0; cc < ncols; cc += ncols_tile)
			for (r = rr; r < rr + nrows_tile; r++)
				for (c = cc; c < cc + ncols_tile ; c++ )
					m_real->rows[r][c] = m->elements[count++];
}



void matrix_print(matrix_t *matrix) 
{
	 unsigned long c,r;

	 for (r = 0; r < matrix->nrows; r++) {
		  for (c = 0; c < matrix->ncols; c++) {
			   printf("%04.4f ", matrix->rows[r][c]);
		  }
		  printf("\n");
	 }
	 return;
}

void vector_print(vector_t *vector) 
{
	 unsigned long i;

	 for (i = 0; i < vector->size; i++) 
		 printf("%lu: %04.4f\n", i, vector->elements[i]);
	 printf("\n");
	 return;
}


/*
 * Print matrix like octave does...
 *
 */ 
void matrix_print_octave(matrix_t *matrix)
{
	 unsigned long c,r;

	 printf("[\n");
	 for (r = 0; ; ) {
		  for (c = 0; ; ) {
			   printf(" %04.4f ", matrix->rows[r][c]);
		   if ( ++c >= matrix->ncols) break;
		   printf(" , ");
		  }
	  if ( ++r >= matrix->nrows ) break;
		  printf(";\n");
	 }
	 printf("\n]\n");
}
#endif
