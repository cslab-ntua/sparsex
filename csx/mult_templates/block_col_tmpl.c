static inline double block_col_${r}x${c}_case(uint8_t **ctl, uint8_t size,
                                              double **values, double **x_curr,
                                              double **y_curr)
{
	register double *y_curr_ = *y_curr;
	register double *x_curr_ = *x_curr;

	for (uint64_t i = 0; i < ${r}; i++) {
		register double yr = 0;
		
		for (uint64_t j = 0; j < ${c}; j++) {
			yr += x_curr_[j] * (**values);
			(*values)++;
		}

		y_curr_[i] += yr;
	}

	return 0;
}
