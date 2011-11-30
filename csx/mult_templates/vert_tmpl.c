static inline double vert${delta}_case(uint8_t **ctl, uint8_t size,
                                        double **values,
                                        double **x_curr, double **y_curr)
{
    register double yr;
    register double xr = **x_curr;
    register uint64_t yi_end = ${delta}*size;

    yr = xr * (**values);
    (*values)++;
    for (uint64_t yi = ${delta}; yi < yi_end; yi += ${delta}) {
        (*y_curr)[yi] += xr * (**values);
        (*values)++;
    }

    return yr;
}
