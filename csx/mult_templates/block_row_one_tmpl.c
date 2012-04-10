static inline double block_row_${r}x${c}_case(uint8_t **ctl, uint8_t size,
                                              double **values, double **x_curr,
                                              double **y_curr)
{
    register double yr = 0;
    register double *x_curr_ = *x_curr;

    for (uint64_t i = 0; i < ${c}; i++) {
        yr += x_curr_[i] * (**values);
        (*values)++;
    }

	return yr;
}
