static inline double delta${bits}_case(uint8_t **ctl, uint8_t size,
                                       double **values, double **x_curr,
                                       double **y_curr)
{
	register double yr;

	${align_ctl}
	yr = (**x_curr) * (**values);
	(*values)++;
	for (uint8_t i = 1; i < size; i++) {
		*x_curr += u${bits}_get(ctl);
		yr += (**x_curr) * (**values);
		(*values)++;
	}

	return yr;
}
