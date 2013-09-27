#include "mat.h"

int main(int argc, char **argv)
{
    libcsx_init();

    /* index_t     rowptr[] = {0,5,6,10,15,18,22,24,29,33,38}; */
    /* index_t     colind[] = {0,1,2,3,8,7,0,1,6,9,0,1,3,5,9,0,1,9,0,1,5,9,2,3,2,3,4,5,7,2,3,4,5,2,3,4,5,9}; */
    /* value_t     values[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26.1,26.2,27,28,29,29.1,29.2,30,31,31.1,31.2,32}; */
    /* index_t     nrows, ncols; */
    /* int         zero_based; */

    /* nrows = 10; ncols = 10; */
    /* zero_based = 1; */
    /* input_t *input = libcsx_mat_create_csr(rowptr, colind, values, */
    /*                                        nrows, ncols, zero_based); */
    input_t *input = libcsx_mat_create_mmf(argv[1]);
//    libcsx_set_tuning_option(input, OPT_REORDER, "ENABLE");
    matrix_t *A = libcsx_mat_tune(input);
    libcsx_mat_destroy_input(input);
    /* libcsx_mat_save(A, NULL); */
    /* libcsx_mat_destroy_tuned(A); */
    /* A = libcsx_mat_restore("csx_file"); */
    vector_t *x = vec_create_random(libcsx_mat_get_ncols(A));
    vector_t *y = vec_create_random(libcsx_mat_get_nrows(A));
    libcsx_matvec_mult(A, 1, x, 1, y);
    printf("\nSpMV output:\n");
    vec_print(y);

    /* Cleanup */
    libcsx_mat_destroy_tuned(A);
    vec_destroy(x);
    vec_destroy(y);
    libcsx_close();

    return 0;
}
