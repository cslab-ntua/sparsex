static inline spx_value_t block_col_${r}x${c}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t **x_curr, spx_value_t **y_curr, spx_value_t scale_f)
{
	register spx_value_t xr = **x_curr;
	register spx_value_t *y_curr_ = *y_curr;

	for (spx_index_t yi = 0; yi < ${r}; yi++) {
		y_curr_[yi] += xr * (**values) * scale_f;
		(*values)++;
	}

	return 0;
}
