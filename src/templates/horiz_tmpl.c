static inline spx_value_t horiz${delta}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t **x_curr, spx_value_t **y_curr, spx_value_t scale_f)
{
    register spx_value_t yr = 0;
    register spx_value_t *x_curr_ = *x_curr;
    register spx_index_t i_end = ${delta}*size;

    for (spx_index_t i = 0; i < i_end; i += ${delta}) {
        yr += x_curr_[i] * (**values);
        (*values)++;
    }

    // Move x_curr to the last element of the pattern
    (*x_curr) += (i_end - ${delta});
	return yr * scale_f;
}
