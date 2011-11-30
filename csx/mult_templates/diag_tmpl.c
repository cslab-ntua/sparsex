static inline double diag${delta}_case(uint8_t **ctl, uint8_t size,
                                       double **values, double **x_curr,
                                       double **y_curr)
{
    register double *y_curr_ = *y_curr;
    register double *x_curr_ = *x_curr;
    register double *end = x_curr_ + ${delta}*size;

    for ( ; x_curr_ < end; x_curr_ += ${delta}, y_curr_ += ${delta}) {
        (*y_curr_) += (*x_curr_) * (**values);
        (*values)++;
    }

    return 0;
}
