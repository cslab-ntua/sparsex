static inline double block_row_${r}x${c}_case(uint8_t **ctl, uint8_t size,
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
	register double curr_x_ = *x_;
	register double curr_ry_ = 0;
	uint64_t i, j;

	for (i = 0; i < ${c}; i++) {
		for (j = 0; j < ${r}; j++) {
			y_[j] += curr_x_ * values_;
			curr_ry_ += rx_[j] * values_;
			(*values)++;
			values_ = **values;
		}
		
		(*ry_) += curr_ry_;
		x_++;
		ry_++;
	    curr_x_ = *x_;
	    curr_ry_ = 0;
	}

	return 0;
}
