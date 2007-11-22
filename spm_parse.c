#include <stdlib.h>
#include <stdio.h>

#include "spm_parse.h"
#include "elem.h"
#include "mmf.h"


void *spm_parse_seq(FILE *mmf, 
                    unsigned long rows_nr, unsigned long cols_nr, 
                    unsigned long nz_nr,
                    spm_parser_t *spm_parser)
{
	parse_state_t parse_st;	
	void *fmt_state, *ret;
	double v;

	parse_st.cur_values_idx = 0;
	parse_st.cur_values = malloc(sizeof(double)*cols_nr);
	if ( !parse_st.cur_values ){
		perror("malloc");
		exit(1);
	}

	if ( !mmf_get_next(mmf, &parse_st.row, &parse_st.col, &v) ){
		fprintf(stderr, "mmf file error: matrix empty\n");
		exit(1);
	}

	spm_parser->initialize(&parse_st, &fmt_state, rows_nr, cols_nr, nz_nr);
	parse_st.cur_values[parse_st.cur_values_idx++] = v;
	parse_st.cur_col = parse_st.cur_unit_col = parse_st.col;
	parse_st.cur_row = parse_st.row;
	
	while ( mmf_get_next(mmf, &parse_st.row, &parse_st.col, &v) ) {
		if ( parse_st.row > parse_st.cur_row ){
			spm_parser->finalize_unit(&parse_st, fmt_state);
			spm_parser->finalize_row(&parse_st, fmt_state);
			parse_st.cur_row = parse_st.row;
			parse_st.cur_unit_col = parse_st.col;
			parse_st.cur_values_idx = 0;
		} else if ( parse_st.col - parse_st.cur_col != 1) {
			spm_parser->finalize_unit(&parse_st, fmt_state);
			parse_st.cur_unit_col = parse_st.col;
			parse_st.cur_values_idx = 0;
		} 
		parse_st.cur_values[parse_st.cur_values_idx++] = v;
		parse_st.cur_col = parse_st.col;
	}
	spm_parser->finalize_unit(&parse_st, fmt_state);
	spm_parser->finalize_row(&parse_st, fmt_state);
	ret = spm_parser->finalize(&parse_st, fmt_state);
	
	free(parse_st.cur_values);

	return ret;
}

#if 0
static void padd_vals_start(parse_state_t *p, unsigned long padd)
{
	unsigned long cur_len = p->cur_values_idx;
	unsigned long i;

	for (i=cur_len-1 ; ; ){
		//printf("i=%lu, i+padd=%lu\n", i, i+padd);
		p->cur_values[i + padd] = p->cur_values[i];
		if ( i == 0){
			break;
		}
		i--;
	}

	for (i=0 ; i<padd ; i++){
		p->cur_values[i] = (elem_t)0;
	}

	p->cur_values_idx += padd;
	p->cur_unit_col -= padd;
}

static void padd_vals_end(parse_state_t *p, unsigned long padd)
{
	unsigned long i;
	unsigned long cur_idx = p->cur_values_idx;
	for (i=0 ; i<padd; i++){
		p->cur_values[cur_idx + i] = (elem_t)0;
	}

	p->cur_values_idx += padd;
}

void vec_prepare_unit(parse_state_t *parse_st, unsigned long align)
{
	unsigned long col_indx = parse_st->cur_unit_col;
	unsigned long padd_s = col_indx & (align-1);
	//printf("col_indx_old:%lu padd_s:%lu\n", col_indx, padd_s);
	if ( padd_s ){
		padd_vals_start(parse_st, padd_s);
	}
	//printf("col_indx_new:%lu\n", parse_st->cur_unit_col);

	unsigned long unit_len = parse_st->cur_values_idx;
	unsigned long padd_e = unit_len & (align-1);
	//printf("len_old:%lu padd_e:%lu\n", unit_len, padd_e);
	if ( padd_e ){
		padd_vals_end(parse_st, padd_e);
	}
	//printf("len_new:%lu\n\n", parse_st->cur_values_idx);
}

/*
 * spm_parse_vec: parsing for vector implementations:
 * Gives units which are:
 *  - size aligned at vec_size
 *  - col idx aligned at vec_size
 *  - uses padding
 */
void *spm_parse_padd(FILE *mmf, 
                    unsigned long rows_nr, unsigned long cols_nr, 
                    unsigned long nz_nr, 
                    spm_parser_t *spm_parser,
                    unsigned long vec_size)
{
	parse_state_t parse_st;	
	void *fmt_state;
	elem_t v;

	parse_st.cur_values_idx = 0;
	parse_st.cur_values = malloc(sizeof(elem_t)*cols_nr);
	if ( !parse_st.cur_values ){
		perror("malloc");
		exit(1);
	}

	if ( !mmf_get_next(mmf, &parse_st.row, &parse_st.col, &v) ){
		fprintf(stderr, "mmf file error: matrix empty\n");
		exit(1);
	}
	spm_parser->initialize(&parse_st, &fmt_state, rows_nr, cols_nr, nz_nr);
	parse_st.cur_values[parse_st.cur_values_idx++] = v;
	parse_st.cur_col = parse_st.cur_unit_col = parse_st.col;
	parse_st.cur_row = parse_st.row;
	
	while ( mmf_get_next(mmf, &parse_st.row, &parse_st.col, &v) ) {
		if ( parse_st.row > parse_st.cur_row ){
			vec_prepare_unit(&parse_st, vec_size);
			spm_parser->finalize_unit(&parse_st, fmt_state);
			spm_parser->finalize_row(&parse_st, fmt_state);
			parse_st.cur_row = parse_st.row;
			parse_st.cur_unit_col = parse_st.col;
			parse_st.cur_values_idx = 0;
		} else {
			unsigned long diff = parse_st.col - parse_st.cur_col;
			if ( diff < vec_size ){
				if ( diff > 1){
					padd_vals_end(&parse_st, diff);
				}
			} else {
				vec_prepare_unit(&parse_st, vec_size);
				spm_parser->finalize_unit(&parse_st, fmt_state);
				parse_st.cur_unit_col = parse_st.col;
				parse_st.cur_values_idx = 0;
			}
		} 
		parse_st.cur_values[parse_st.cur_values_idx++] = v;
		parse_st.cur_col = parse_st.col;
	}

	vec_prepare_unit(&parse_st, vec_size);
	spm_parser->finalize_unit(&parse_st, fmt_state);
	spm_parser->finalize_row(&parse_st, fmt_state);
	spm_parser->finalize(&parse_st, fmt_state);
}
#endif
