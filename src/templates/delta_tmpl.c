static inline spx_value_t delta${bits}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t **x_curr, spx_value_t **y_curr, spx_value_t scale_f)
{
	register spx_value_t yr;

	${align_ctl}
	yr = (**x_curr) * (**values);
	(*values)++;
	for (uint8_t i = 1; i < size; i++) {
		*x_curr += u${bits}_get(ctl);
		yr += (**x_curr) * (**values);
		(*values)++;
	}

	return yr * scale_f;
}
