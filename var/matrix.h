#ifndef __MATRIX_H__
#define __MATRIX_H__

#if 0
#ifndef ELEM_SIZE
#define ELEM_SIZE 8
#endif

#if ELEM_SIZE == 4
#define ELEM_TYPE float
#define ELEM_TYPE_FLOAT
#elif ELEM_SIZE == 8
#define ELEM_TYPE double
#define ELEM_TYPE_DOUBLE
#else
#error "Invalid ELEM_SIZE"
#endif

typedef ELEM_TYPE elem_t;
#endif

enum {
	ELEM_FLOAT = 1,
	ELEM_DOUBLE
};

static inline int elem_size(unsigned type)
{
	switch (type){
		case ELEM_FLOAT:
		return sizeof(float);

		case ELEM_DOUBLE:
		return sizeof(double);
		
		default:
		fprintf(stderr, "%s: uknown type: %u\n", __FUNCTION__, type);
		exit(1);
	}
}


typedef struct matrix_s {
	void *elements;
	void **rows;
	unsigned long nrows, ncols;
	unsigned elem_type;
} matrix_t;

typedef struct vector_s {
	void *elements;
	unsigned long size;
	unsigned elem_type;
} vector_t;

matrix_t *matrix_create(unsigned long nrows, unsigned long ncols, unsigned elem_type);
void matrix_destroy(matrix_t *matrix);
extern vector_t *vector_create(unsigned long size, unsigned elem_type);
extern void vector_destroy(vector_t *vector);
#ifdef NUMA
extern vector_t *numa_vector_create(unsigned long size, int memnode_number);
extern void numa_vector_destroy(vector_t *vector);
#endif
/*
 *	Initialization routines
 */
void matrix_init(matrix_t *matrix, double val);
void matrix_init_unit(matrix_t *matrix);
void matrix_init_inc(matrix_t *matrix, int max_val);
//extern void matrix_init_rand(matrix_t *matrix, int max_val);
void matrix_init_rand_range(matrix_t *matrix, double max, double min);

void vector_init(vector_t *vector, double val);
//void vector_init_rand(vector_t *vector, int max_val);
//void vector_copy(vector_t *src,vector_t *dest);
//extern void vector_split2(vector_t *vector, vector_t *first, vector_t *second, unsigned long where);
void vector_init_rand_range(vector_t *vector, double max, double min);

#if 0
/*
 * Layout routines
 */ 
extern void matrix_row2block_ZZ(matrix_t *m_real, matrix_t *m, unsigned int nrows_tile, unsigned int ncols_tile);
extern void matrix_block2row_ZZ(matrix_t *m, matrix_t *m_real, unsigned int nrows_tile, unsigned int ncols_tile);

/*
 * Utility routines
 */ 
extern void matrix_print(matrix_t *matrix);
extern void matrix_print_octave(matrix_t *matrix);

extern void vector_print(vector_t *vector);
#endif

#endif
