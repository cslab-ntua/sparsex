static inline spx_value_t block_row_${r}x${c}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t **x_curr, spx_value_t **y_curr, spx_value_t scale_f)
{
    register spx_value_t yr = 0;
    register spx_value_t *x_curr_ = *x_curr;

    for (spx_index_t i = 0; i < ${c}; i++) {
        yr += x_curr_[i] * (**values);
        (*values)++;
    }

	return yr * scale_f;
}
