static inline double rdiag${delta}_case(uint8_t **ctl, uint8_t size,
                                        double **values, double **x_curr,
                                        double **y_curr)
{
	register double *x_curr_ = *x_curr;
	register double *y_curr_ = *y_curr;
    register uint64_t i_end = ${delta}*size;

    for (uint64_t i = 0; i < i_end; i += ${delta}) {
        y_curr_[i] += (*x_curr_) * (**values);
        x_curr_ -= ${delta};
        (*values)++;
    }

	return 0;
}
