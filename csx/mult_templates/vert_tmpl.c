static inline double vert${delta}_case(uint8_t **ctl, uint8_t size,
                                       double **values, double **x_curr,
                                       double **y_curr)
{
	register double xr = **x_curr;
	register uint64_t yi_end = ${delta}*size;
	register double *y_curr_ = *y_curr;

	for (uint64_t yi = 0; yi < yi_end; yi += ${delta}) {
		y_curr_[yi] += xr * (**values);
		(*values)++;
	}

	return 0;
}
