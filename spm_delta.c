#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "spm_parse.h"
#include "spm_delta.h"
#include "beanbag.h"
#include "mmf.h"
#include "method.h"

#if 0

#define MIN(a,b) ( (a < b) ? a : b )
#define MAX(a,b) ( (a > b) ? a : b )

#ifndef ELEM_TYPE
#warning "ELEM_TYPE undefined, defining to double"
#define ELEM_TYPE double
#endif

#define SPM_DELTA_CON3_(a,b,c) a ## b ## c
#define SPM_DELTA_CON3(a,b,c) SPM_DELTA_CON3_(a,b,c)
#define SPM_DELTA_NAME(name) SPM_DELTA_CON3(spm_delta_, ELEM_TYPE, name)

static inline int is_aligned(unsigned long ptr, unsigned int align)
{
	return !(ptr & (align-1) );
}

void spm_delta_ctl_print(unsigned char flags, unsigned int size)
{
	printf("nr:%d sp:%d nr_seq:%d ci_size:%d size:%d\n",
               !!spm_delta_fl_isnr(flags),
               !!spm_delta_fl_issp(flags),
               !!spm_delta_fl_isnrseq(flags),
	       spm_delta_fl_getcisize(flags),
	       size);
}

spm_delta_t *spm_delta_create(unsigned long rows_nr, unsigned long cols_nr,
                              unsigned long nz_nr)
{
	spm_delta_t *delta;

	delta = malloc(sizeof(spm_delta_t));
	if ( !delta ){
		perror("spm_delta_create: malloc");
		exit(1);
	}

	delta->values = malloc(sizeof(elem_t)*nz_nr);
	if ( !delta->values ){
		perror("spm_delta_create: malloc");
		exit(1);
	}

	delta->nnz = nz_nr;
	delta->ncols = cols_nr;
	delta->nrows = rows_nr;

	return delta;
}

struct delta_sp_state_s {
	elem_t          *values;
	unsigned long   *ci_deltas;
	unsigned char   ci_max;
	unsigned long   idx;
	unsigned long   idx_consumed;
};
typedef struct delta_sp_state_s delta_sp_state_t;

struct delta_state_s {
	spm_delta_t      *delta;
	beanbag_t        *ctl_bb;
	unsigned long    val_indx, col_indx, row_indx, curval_indx;
	char             new_row;

	delta_sp_state_t sp;
};
typedef struct delta_state_s delta_state_t;

static inline unsigned long delta_colind_get_delta(delta_state_t *delta_state, unsigned long colind)
{
	unsigned long ret;
	if ( colind < delta_state->col_indx){
		fprintf(stderr, "delta_colind_get_delta: %lu-%lu=%lu\n", colind, delta_state->col_indx, colind-delta_state->col_indx);
		exit(1);
	}
	ret = colind-delta_state->col_indx;
	delta_state->col_indx = colind;

	return ret;
}


static void delta_init(parse_state_t *parse_state, void **delta_state_ptr,
                       unsigned long rows_nr, unsigned long cols_nr,
                       unsigned long nz_nr)
{
	spm_delta_t   *delta;
	delta_state_t *delta_state;

	delta = spm_delta_create(rows_nr, cols_nr, nz_nr);
	delta_state = malloc(sizeof(delta_state_t));
	if ( !delta_state ){
		perror("delta_init: malloc");
		exit(1);
	}

	delta_state->sp.values = malloc(sizeof(elem_t)*cols_nr);
	delta_state->sp.ci_deltas = malloc(sizeof(unsigned long)*cols_nr);
	if ( !delta_state->sp.values || !delta_state->sp.ci_deltas){
		perror("delta_init: mallloc");
		exit(1);
	}

	delta_state->delta = delta;
	delta_state->ctl_bb = beanbag_create(sizeof(char), 4096);
	delta_state->val_indx = 0;
	delta_state->col_indx = 0;
	delta_state->row_indx = 0;
	delta_state->new_row = 0;

	delta_state->sp.idx = 0;
	delta_state->sp.idx_consumed = 0;
	delta_state->sp.ci_max = SPM_DELTA_CISIZE_U8;

	*delta_state_ptr = delta_state;
}

#define LONGUC_SHIFT (7)
static void delta_bb_append_ul(beanbag_t *bb, unsigned long val)
{
	unsigned char *uc;
	const unsigned shift = LONGUC_SHIFT;

	for (;;) {
		uc = beanbag_alloc(bb);
		*uc = (val & ((1<<shift) - 1));
		if ( val < (1<<shift) ){
			break;
		}
		*uc |= (1<<shift);
		val >>= shift;
	}
}
#if 0
static unsigned long delta_uc_to_ul(unsigned char **uc_ptr)
{
	unsigned sh=0;
	unsigned long val = 0, _val;
	const unsigned shift = LONGUC_SHIFT;

	for (;;){
		_val = **uc_ptr;
		(*uc_ptr)++;
		val |= ( (_val & (~(1<<shift))) << sh);
		if ( !(_val & (1<<shift)) ) {
			break;
		}
		sh += shift;
	}
	return val;
}
#endif


#if 0
#define LONGUC_SHIFT (7)
#define LONGUS_SHIFT (15)
static void delta_bb_append_ul(beanbag_t *bb, unsigned long val)
{
	unsigned short *us;
	const unsigned shift = LONGUS_SHIFT;

	for (;;) {
		us = beanbag_alloc_nr(bb, 2);
		*us = (val & ((1<<shift) - 1));
		if ( val < (1<<shift) ){
			break;
		}
		*us |= (1<<shift);
		val >>= shift;
	}
}

static unsigned long delta_uc_to_ul(unsigned char **uc_ptr)
{
	unsigned sh=0;
	register unsigned long val = 0, _val;
	const unsigned shift = LONGUS_SHIFT;

	_val = *((unsigned short *)(*uc_ptr));
	if ( _val < 16384) {
		*uc_ptr += 2;
		return _val;
	}

	val = _val - 128;
	_val = *(*uc_ptr + 1);
	val += (_val<<shift);
	if ( _val < 128 ){
		*uc_ptr += 2;
		return val;
	}

	return val;
}
#endif

static void delta_do_new_row(parse_state_t *parse_state,
                             delta_state_t *delta_state,
			     unsigned char *ctl_flags)
{
	unsigned long delta_row = parse_state->cur_row - delta_state->row_indx;
	spm_delta_fl_setnr(ctl_flags);
	if ( delta_row == 1){
		spm_delta_fl_setnrseq(ctl_flags);
	} else {
		delta_bb_append_ul(delta_state->ctl_bb, delta_row);
	}
	delta_state->new_row = 0;
	delta_state->row_indx = parse_state->cur_row;
}


static void delta_do_finalize_unit(parse_state_t *parse_state, delta_state_t *delta_state,
                                   unsigned char size)
{
	unsigned char *ctl_flags, *ctl_size, *ctl;
	elem_t *v, *v_cur;
	unsigned long i, delta_col;
	spm_delta_t *delta = delta_state->delta;

	/* allocate and create proper ctl */
	ctl = beanbag_alloc_nr(delta_state->ctl_bb, 2);
	ctl_flags = &ctl[0];
	ctl_size = &ctl[1];

	*ctl_flags = 0;
	*ctl_size = size;

	if ( delta_state->new_row ){
		delta_do_new_row(parse_state, delta_state, ctl_flags);
	}

	//printf("DE: size:%u rind:%lu ", size, delta_state->row_indx);
	//spm_delta_ctl_print(*ctl_flags, *ctl_size);

	/* add column index */
	delta_col = delta_colind_get_delta(delta_state, parse_state->cur_unit_col);
	delta_bb_append_ul(delta_state->ctl_bb, delta_col);

	/* add values */
	v = &delta->values[delta_state->val_indx];
	v_cur = &parse_state->cur_values[delta_state->curval_indx];
	for (i=0 ; i < size ; i++){
		v[i] = v_cur[i];
	}
	parse_state->cur_unit_col += size;
	delta_state->val_indx += size;
	delta_state->col_indx += size;
}

static void delta_do_sp_finalize(parse_state_t *parse_state,
                                 delta_state_t *delta_state,
                                 delta_sp_state_t *sp_state,
                                 unsigned char size)
{
	spm_delta_t *delta = delta_state->delta;
	unsigned char *ctl, *ctl_flags, *ctl_size, ci_deltas_size;
	void *ci_deltas;
	elem_t *v_from, *v_to;
	unsigned long i, sp_idx, *ci_deltas_from;

	ci_deltas_size = spm_delta_cisize_size(sp_state->ci_max);
	ctl = beanbag_alloc_nr(delta_state->ctl_bb, 2);
	ctl_flags = &ctl[0];
	ctl_size = &ctl[1];


	*ctl_flags = 0;
	spm_delta_fl_setsp(ctl_flags);
	spm_delta_fl_setcisize(ctl_flags, sp_state->ci_max);
	*ctl_size = size;

	if ( delta_state->new_row ){
		delta_do_new_row(parse_state, delta_state, ctl_flags);
	}

	//printf("SP: size:%u rind:%lu ci_max:%u ", size, delta_state->row_indx, sp_state->ci_max);
	//spm_delta_ctl_print(*ctl_flags,*ctl_size);

	sp_idx = sp_state->idx_consumed;

	ci_deltas = beanbag_alloc_nr(delta_state->ctl_bb,ci_deltas_size*size);
	v_from = &sp_state->values[sp_idx];
	v_to = &delta->values[delta_state->val_indx];
	ci_deltas_from = &sp_state->ci_deltas[sp_idx];
	for (i=0; i<size; i++){
		v_to[i] = v_from[i];
	}
	spm_delta_cisize_copy(ci_deltas, ci_deltas_from, size, sp_state->ci_max);
	delta_state->val_indx += size;
	sp_state->idx_consumed += size;
}

static void delta_sp_finalize(parse_state_t *parse_state,
                              delta_state_t *delta_state)
{
	delta_sp_state_t *sp_state = &delta_state->sp;
	unsigned char size;
	unsigned long len = sp_state->idx;

	sp_state->idx_consumed = 0;
	while ( len ){
		size = MIN(SPM_DELTA_CTL_SIZE_MAX, len);
		delta_do_sp_finalize(parse_state, delta_state, sp_state, size);
		len -= size;
	}
	sp_state->idx = 0;
	sp_state->ci_max = SPM_DELTA_CISIZE_U8;
}

static void delta_sp_add_element(delta_state_t *delta_state, unsigned long colind, elem_t val)
{
	delta_sp_state_t *sp_state = &delta_state->sp;
	unsigned long idx = sp_state->idx;
	unsigned long ci_delta = delta_colind_get_delta(delta_state, colind);
	unsigned char ci_delta_size = spm_delta_cisize(ci_delta);

	//printf("%s: colind=%lu val=%lf\n", __FUNCTION__, colind, val);
	if (ci_delta_size > sp_state->ci_max){
		sp_state->ci_max = ci_delta_size;
	}

	sp_state->values[idx] = val;
	sp_state->ci_deltas[idx] = ci_delta;
	sp_state->idx++;
}

static void delta_sp_add_unit(parse_state_t *parse_state,
                              delta_state_t *delta_state,
                              unsigned long unit_len)
{
	unsigned long i;
	unsigned long col_indx = parse_state->cur_unit_col;

	for ( i=0; i < unit_len; i++){
		delta_sp_add_element(delta_state, col_indx+i, parse_state->cur_values[i]);
	}
}

#define UNIT_MIN_SEQ_LEN 4
static void delta_finalize_unit(parse_state_t *parse_state, void *delta_st)
{
	delta_state_t  *delta_state = delta_st;
	unsigned long  unit_len = parse_state->cur_values_idx;

	assert(unit_len);
	if ( unit_len > UNIT_MIN_SEQ_LEN ) {
		if ( delta_state->sp.idx ){
			/* there is a sparse unit waiting to be finalized */
			delta_sp_finalize(parse_state, delta_state);
		}
		/* create the appropriate number of dense units */
		delta_state->curval_indx = 0;
		while ( unit_len ){
			unsigned char size = MIN(SPM_DELTA_CTL_SIZE_MAX, unit_len);
			delta_do_finalize_unit(parse_state, delta_state, size);
			unit_len -= size;
			delta_state->curval_indx += size;
		}
	} else {
		/* add data to the current sparse unit */
		delta_sp_add_unit(parse_state, delta_state, unit_len);
	}
}

static void delta_finalize_row(parse_state_t *parse_state, void *delta_st)
{
	delta_state_t  *delta_state = delta_st;
	if ( delta_state->sp.idx){
		/* there is a sparse unit waiting to be finalized */
		delta_sp_finalize(parse_state, delta_state);
	}
	delta_state->col_indx = 0;
	/* set the new row flag */
	delta_state->new_row = 1;
}

static void *delta_final(parse_state_t *parse_state, void *delta_st)
{
	delta_state_t  *delta_state = delta_st;
	spm_delta_t    *delta = delta_state->delta;

	if ( delta_state->sp.idx ){
		delta_sp_finalize(parse_state, delta_state);
	}

	delta->nnz = delta_state->val_indx;
	printf("ctl size: %lu\n", beanbag_size(delta_state->ctl_bb));
	delta->ctl = beanbag_destroy(delta_state->ctl_bb);
	free(delta_state->sp.values);
	free(delta_state->sp.ci_deltas);
	free(delta_state);

	return delta;
}

static spm_parser_t delta_parser = {
	.initialize = delta_init,
	.finalize_unit = delta_finalize_unit,
	.finalize_row = delta_finalize_row,
	.finalize = delta_final,
};

spm_delta_t *spm_delta_init_mmf(char *mmf_file,
                                unsigned long *rows_nr, unsigned long *cols_nr,
                                unsigned long *nz_nr)
{
	int ret;
	spm_delta_t *delta;
	FILE *f;

	f = mmf_init(mmf_file, rows_nr, cols_nr, nz_nr);
	delta = spm_parse_seq(f, *rows_nr, *cols_nr, *nz_nr, &delta_parser);
	fclose(f);

	return delta;
}

#endif


#if 0
void spm_delta_multiply(spm_delta_t *delta, vector_t *in, vector_t *out)
{
	elem_t *x = in->elements, *y = out->elements, *values = delta->values;
	register unsigned long y_indx=0, v_indx=0;
	register elem_t yr = 0;
	register elem_t *v = values;
	register const elem_t *v_end = values + delta->nnz;
	register elem_t *myx = x;
	register unsigned char *uc = delta->ctl;

	for (;;) {
		register unsigned char flags = uc_get(uc);
		register unsigned char size = uc_get(uc);
		//spm_delta_ctl_print(flags, size);
		if ( spm_delta_fl_isnr(flags) ){
			y[y_indx] = yr;
			if ( spm_delta_fl_isnrseq(flags) ){
				y_indx++;
			} else {
				y_indx += uc_get_ul(uc);
			}
			//printf("--new_row: (result for %lu: %lf)\n", y_indx, yr);
			yr = 0;
			myx = x;
		}
		//printf("y_indx: %lu x_indx: %lu\n", y_indx, (unsigned long)(myx-x));

		switch ( flags & SPM_DELTA_FL_UNIT_MASK ){
			register unsigned long i;

			case SPM_DELTA_FL_UNIT_DENSE:
			__asm__ __volatile__ (" # DE loop ");
			myx += uc_get_ul(uc);
			for ( i = 0; i<size; i++){
				yr += (*v) * (*myx);
				//printf("\t[DE] y_indx:%lu vindx=%lu xindx=%lu (%lf %lf %lf)\n", y_indx, (unsigned long)(v - values), (unsigned long)(myx - x), *v, *myx, yr);
				v++;
				myx++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_UC:
			for ( i = 0; i<size; i++){
				myx += uc_get(uc);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_US:
			for ( i = 0; i<size; i++){
				myx += us_get(uc);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_UI:
			for ( i = 0; i<size; i++){
				myx += ui_get(uc);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			case SPM_DELTA_FL_UNIT_SP_UL:
			for ( i = 0; i<size; i++){
				myx += ul_get(uc);
				yr += (*v) * (*myx);
				v++;
			}
			break;

			default:
			printf("Uknown flags: %u unit=%u\n", flags, flags & SPM_DELTA_FL_UNIT_MASK);
			assert(0);
		}

		if ( v >= v_end ){
			break;
		}
	}

	y[y_indx] = yr;
	return;
}
#endif
//METHOD_INIT(spm_delta_multiply, spm_delta_init_mmf)
