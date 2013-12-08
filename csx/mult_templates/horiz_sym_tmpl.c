static inline double horiz${delta}_case(uint8_t **ctl, uint8_t size,
                                        double **values, double *x, double *y,
                                        double *cur, uint64_t *x_indx,
                                        uint64_t *y_indx, double scale_f)
{
	register double yr = 0;
	register double values_ = **values;
	register uint64_t x_indx_ = *x_indx;
	register double *x_ = x + x_indx_;
	register double *ry_ = cur + x_indx_;
	register uint64_t y_indx_ = *y_indx;
	register double rx_ = x[y_indx_];
	register uint64_t i_end = ${delta}*size;

	for (uint64_t i = 0; i < i_end; i += ${delta}) {
		yr += x_[i] * values_;
		ry_[i] += rx_ * values_ * scale_f;
		(*values)++;
		values_ = **values;
	}

	// Move x_indx to the last element of the pattern
	(*x_indx) += i_end - ${delta};
	return yr * scale_f;
}
