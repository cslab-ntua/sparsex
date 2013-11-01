static inline double block_col_${r}x${c}_case(uint8_t **ctl, uint8_t size,
                                              double **values, double **x_curr,
                                              double **y_curr, double scale_f)
{
	register double xr = **x_curr;
	register double *y_curr_ = *y_curr;

	for (uint64_t yi = 0; yi < ${r}; yi++) {
		y_curr_[yi] += xr * (**values) * scale_f;
		(*values)++;
	}

	return 0;
}
