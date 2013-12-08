#include "libcsx.h"

int main(int argc, char **argv)
{
    /* index_t     rowptr[] = {0,5,6,10,15,18,22,24,29,33,38}; */
    /* index_t     colind[] = {0,1,2,3,8,7,0,1,6,9,0,1,3,5,9,0,1,9,0,1,5,9,2,3,2,3,4,5,7,2,3,4,5,2,3,4,5,9}; */
    /* value_t     values[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26.1,26.2,27,28,29,29.1,29.2,30,31,31.1,31.2,32}; */
    /* index_t     nrows, ncols; */

    /* nrows = 10; ncols = 10; */
    /* input_t *input = libcsx_mat_create_csr(rowptr, colind, values, */
    /*                                        nrows, ncols, 
                                              INDEXING_ZERO_BASED); */

    libcsx_init();
    printf("Loading matrix from file...\n");
    input_t *input = libcsx_mat_create_mmf(argv[1]);

    printf("Transforming to CSX...\n");
    libcsx_set_options_from_env();
    /* libcsx_set_option("libcsx.matrix.symmetric", "true"); */
    matrix_t *A = libcsx_mat_tune(input);
    /* printf("Dumping to binary...\n"); */
    /* libcsx_mat_save(A, NULL); */
    /* libcsx_mat_destroy_tuned(A); */
    /* A = libcsx_mat_restore("csx_file"); */

    printf("Vector creation...\n");
    vector_t *x = libcsx_vec_create_random(libcsx_mat_get_ncols(A), A);
    vector_t *y = libcsx_vec_create_random(libcsx_mat_get_nrows(A), A);

    printf("Running SpMV kernel...\n");
    libcsx_matvec_mult(A, 0.58, x, 0.1, y);

    /* printf("SpMV output:\n"); */
    /* libcsx_vec_print(y); */

    /* libcsx_init(); */
    /* printf("Loading matrix from file...\n"); */
    /* input_t *input = libcsx_mat_create_mmf(argv[1]); */
    /* printf("Transforming to CSX...\n"); */
    /* libcsx_set_options_from_env(); */
    /* /\* libcsx_set_option("libcsx.matrix.symmetric", "true"); *\/ */
    /* matrix_t *A = libcsx_mat_tune(input, OP_REORDER); */

    /* /\* printf("Dumping to binary...\n"); *\/ */
    /* /\* libcsx_mat_save(A, NULL); *\/ */
    /* /\* libcsx_mat_destroy_tuned(A); *\/ */
    /* /\* A = libcsx_mat_restore("csx_file"); *\/ */

    /* printf("Vector creation...\n"); */
    /* vector_t *x = libcsx_vec_create_random(libcsx_mat_get_ncols(A), A); */
    /* vector_t *y = libcsx_vec_create_random(libcsx_mat_get_nrows(A), A); */
    /* perm_t *p = libcsx_mat_get_perm(A); */
    /* vector_t *permuted_x = libcsx_vec_reorder(x, p); */
    /* vector_t *permuted_y = libcsx_vec_reorder(y, p); */

    /* printf("Running SpMV kernel...\n"); */
    /* libcsx_matvec_mult(A, 0.1, permuted_x, 0.41, permuted_y); */

    /* printf("SpMV output:\n"); */
    /* vector_t *tmp_y = libcsx_vec_inv_reorder(permuted_y, p); */
    /* vec_copy(tmp_y, y); */
    /* vec_destroy(tmp_y); */
    /* libcsx_vec_print(y); */
    /* libcsx_vec_destroy(permuted_x); */
    /* libcsx_vec_destroy(permuted_y); */

    /* Cleanup */
    libcsx_mat_destroy_input(input);
    libcsx_mat_destroy_tuned(A);
    libcsx_vec_destroy(x);
    libcsx_vec_destroy(y);

    return 0;
}
