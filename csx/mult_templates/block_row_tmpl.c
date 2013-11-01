static inline double block_row_${r}x${c}_case(uint8_t **ctl, uint8_t size,
                                              double **values, double **x_curr,
                                              double **y_curr, double scale_f)
{
	register double *y_curr_ = *y_curr;
	register double *x_curr_ = *x_curr;

	for (uint64_t i = 0; i < ${c}; i++) {
		register double xr = x_curr_[i];
		
		for (uint64_t j = 0; j < ${r}; j++) {
			y_curr_[j] += xr * (**values) * scale_f;
			(*values)++;
		}
	}

	return 0;
}
