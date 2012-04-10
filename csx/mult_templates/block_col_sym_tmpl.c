static inline double block_col_${r}x${c}_case(uint8_t **ctl, uint8_t size,
                                              double **values, double *x,
                                              double *y, double *cur,
                                              uint64_t *x_indx,
                                              uint64_t *y_indx)
{
	register double values_ = **values;
	register uint64_t x_indx_ = *x_indx;
	register double *x_ = x + x_indx_;
	register double *ry_ = cur + x_indx_;
	register uint64_t y_indx_ = *y_indx;
	register double *y_ = y + y_indx_;
	register double *rx_ = x + y_indx_;
	register double curr_y_ = 0;
	register double curr_rx_ = *rx_;
	uint64_t i, j;

	for (i = 0; i < ${r}; i++) {
		for (j = 0; j < ${c}; j++) {
			curr_y_ += x_[j] * values_;
			ry_[j] += curr_rx_ * values_;
			(*values)++;
			values_ = **values;
		}
		
		(*y_) += curr_y_;
		y_++;
		rx_++;
	    curr_y_ = 0;
	    curr_rx_ = *rx_;
	}

	return 0;
}
