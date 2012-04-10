static inline double vert${delta}_case(uint8_t **ctl, uint8_t size,
                                       double **values, double *x, double *y,
                                       double *cur, uint64_t *x_indx,
                                       uint64_t *y_indx)
{
	register double values_ = **values;
	register uint64_t x_indx_ = *x_indx;
	register double x_ = x[x_indx_];
	register double ry_ = 0;
	register uint64_t y_indx_ = *y_indx;
	register double *y_ = y + y_indx_;
	register double *rx_ = x + y_indx_;
	register uint64_t i_end = ${delta} * size;

	for (uint64_t i = 0; i < i_end; i += ${delta}) {
		y_[i] += x_ * values_;
		ry_ += rx_[i] * values_;
		(*values)++;
		values_ = **values;
	}
    cur[x_indx_] += ry_;
    
    return 0;
}
