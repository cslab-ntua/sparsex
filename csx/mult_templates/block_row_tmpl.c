static inline double block_row_${r}x${c}_case(uint8_t **ctl, uint8_t size,
                                              double **values,
                                              double **x_curr, double **y_curr)
{
    register double *y_curr_ = *y_curr;
    register double *x_curr_ = *x_curr;

    for (uint8_t i = 0; i < ${c}; i++) {
        register double xr = x_curr_[i];
        ${intern_loop}
        for (uint8_t j = 0; j < ${r}; j++) {
            y_curr_[j] += xr * (**values);
            (*values)++;
        }
    }
    
    return 0;
}
