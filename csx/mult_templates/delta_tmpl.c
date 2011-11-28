static inline double delta${bits}_case(uint8_t **ctl, uint8_t size,
                                       double **values,
                                       double **x_curr, double **y_curr)
{
    register double yr = 0;
    for (uint8_t i = 0; i < size; i++) {
        *x_curr += u${bits}_get(ctl);
//        ctl_print(*ctl, 0, 1, "inside delta");
        yr += (**x_curr) * (**values);
        (*values)++;
    }

    return yr;
}
