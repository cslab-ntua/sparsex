static inline double delta${bits}_case(uint8_t **ctl, uint8_t size,
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
    register uint64_t jmp;

    ${align_ctl}
	yr += (*x_) * values_;
	(*ry_) += rx_ * values_ * scale_f;
	(*values)++;
	for (uint8_t i = 1; i < size; i++) {
		jmp = u${bits}_get(ctl);
		x_ += jmp;
		ry_ += jmp;
		(*x_indx) += jmp;
		values_ = **values;
		yr += (*x_) * values_;
		(*ry_) += rx_ * values_ * scale_f;
		(*values)++;
	}

	return yr * scale_f;
}
