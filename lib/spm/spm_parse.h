#ifndef __SPM_PARSE_H__
#define __SPM_PARSE_H__

#include <stdio.h>

struct parse_state_s {
	double *cur_values;
	unsigned long cur_values_idx;
	unsigned long cur_row, cur_col, cur_unit_col;
	unsigned long row, col;
};
typedef struct parse_state_s parse_state_t;

typedef void parser_fn_t(parse_state_t *pstate, void *own_state);
typedef void parser_init_fn_t(parse_state_t *pstate, void **own_state, unsigned long rows_nr, unsigned long cols_nr, unsigned long nz_nr);
typedef void *parser_fin_fn_t(parse_state_t *pstate, void *own_state);

struct spm_parser_s {
	parser_init_fn_t   *initialize;
	parser_fin_fn_t    *finalize;
	parser_fn_t        *finalize_unit, *finalize_row;
};
typedef struct spm_parser_s spm_parser_t;

void *spm_parse_seq(FILE *mmf,
                    unsigned long rows_nr, unsigned long cols_nr,
                    unsigned long nz_nr,
                    spm_parser_t *parser);

#endif /* __SPM_PARSE_H__ */
