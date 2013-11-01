static inline double horiz${delta}_case(uint8_t **ctl, uint8_t size,
                                        double **values, double **x_curr,
                                        double **y_curr, double scale_f)
{
    register double yr = 0;
    register double *x_curr_ = *x_curr;
    register uint64_t i_end = ${delta}*size;

    for (uint64_t i = 0; i < i_end; i += ${delta}) {
        yr += x_curr_[i] * (**values);
        (*values)++;
    }

    // Move x_curr to the last element of the pattern
    (*x_curr) += (i_end - ${delta});
	return yr * scale_f;
}
