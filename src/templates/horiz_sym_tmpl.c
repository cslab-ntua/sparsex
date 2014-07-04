static inline spx_value_t horiz${delta}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t *x, spx_value_t *y, spx_value_t *cur,
    spx_index_t *x_indx, spx_index_t *y_indx, spx_value_t scale_f)
{
	register spx_value_t yr = 0;
	register spx_value_t values_ = **values;
	register spx_index_t x_indx_ = *x_indx;
	register spx_value_t *x_ = x + x_indx_;
	register spx_value_t *ry_ = cur + x_indx_;
	register spx_index_t y_indx_ = *y_indx;
	register spx_value_t rx_ = x[y_indx_];
	register spx_index_t i_end = ${delta}*size;

	for (spx_index_t i = 0; i < i_end; i += ${delta}) {
		yr += x_[i] * values_;
		ry_[i] += rx_ * values_ * scale_f;
		(*values)++;
		values_ = **values;
	}

	// Move x_indx to the last element of the pattern
	(*x_indx) += i_end - ${delta};
	return yr * scale_f;
}
