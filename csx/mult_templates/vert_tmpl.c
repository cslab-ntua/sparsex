static inline double vert${delta}_case(uint8_t **ctl, uint8_t size,
                                        double **values,
                                        double **x_curr, double **y_curr)
{
    register double yr;

    yr = (**x_curr) * (**values);
    (*values)++;
    for (uint64_t yi = ${delta}; yi < ${delta}*size; yi += ${delta}) {
        (*y_curr)[yi] += (**x_curr) * (**values);
        (*values)++;
    }

    return yr;
}
