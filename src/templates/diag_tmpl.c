static inline spx_value_t diag${delta}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t **x_curr, spx_value_t **y_curr, spx_value_t scale_f)
{
	register spx_value_t *y_curr_ = *y_curr;
	register spx_value_t *x_curr_ = *x_curr;
	register spx_index_t i_end = ${delta}*size;

	for (spx_index_t i = 0; i < i_end; i += ${delta}) {
		y_curr_[i] += x_curr_[i] * (**values) * scale_f;
		(*values)++;
	}

	return 0;
}
