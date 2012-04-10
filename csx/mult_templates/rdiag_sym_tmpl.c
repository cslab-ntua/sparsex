static inline double rdiag${delta}_case(uint8_t **ctl, uint8_t size,
                                        double **values, double *x, double *y,
                                        double *cur, uint64_t *x_indx,
                                        uint64_t *y_indx)
{
	register double values_ = **values;
	register uint64_t i_end = ${delta} * size;
	register uint64_t x_indx_ = *x_indx;
	register double *x_ = x + x_indx_ - i_end;
	register double *ry_ = cur + x_indx_ - i_end;
	register uint64_t y_indx_ = *y_indx;
	register double *y_ = y + y_indx_;
	register double *rx_ = x + y_indx_;
	uint64_t i, j;

	for (i = 0, j = i_end; i < i_end; i += ${delta}, j -= ${delta}) {
		y_[i] += x_[j] * values_;
		ry_[j] += rx_[i] * values_;
		(*values)++;
		values_ = **values;
	}

	return 0;
}
