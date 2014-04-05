#include "sparsex.h"

int main(int argc, char **argv)
{
    /* spx_index_t     rowptr[] = {0,5,6,10,15,18,22,24,29,33,38}; */
    /* spx_index_t     colind[] = {0,1,2,3,8,7,0,1,6,9,0,1,3,5,9,0,1,9,0,1,5,9,2,3,2,3,4,5,7,2,3,4,5,2,3,4,5,9}; */
    /* spx_value_t     values[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26.1,26.2,27,28,29,29.1,29.2,30,31,31.1,31.2,32}; */
    /* spx_index_t nrows, ncols; */

    /* nrows = 10; ncols = 10; */
    /* input_t *input = spx_mat_create_csr(rowptr, colind, values, */
    /*                                        nrows, ncols, 
                                              INDEXING_ZERO_BASED); */

    spx_init();

    printf("Loading matrix from file...\n");
    spx_input_t *input = spx_input_load_mmf(argv[1]);

    printf("Transforming to CSX...\n");
    spx_options_set_from_env();

    /* spx_option_set("spx.matrix.symmetric", "true"); */
    spx_matrix_t *A = spx_mat_tune(input);

    printf("Creating x and y vectors...\n");
    spx_partition_t *parts = spx_mat_get_parts(A);
    spx_vector_t *x = spx_vec_create_random(spx_mat_get_ncols(A), parts);
    spx_vector_t *y = spx_vec_create_random(spx_mat_get_nrows(A), parts);
    spx_vector_t *perm_x = spx_vec_create(spx_mat_get_ncols(A), parts);
    spx_vector_t *perm_y = spx_vec_create(spx_mat_get_nrows(A), parts);
    vec_copy(x, perm_x);
    vec_copy(y, perm_y);

    printf("Running the SpMV kernel...\n");
    spx_matvec_mult(0.1, A, x, 0.41, y);

    printf("Result: ");
    /* spx_vec_print(y); */

    printf("Transforming to CSX with reordering enabled...\n");
    spx_matrix_t *B = spx_mat_tune(input, OP_REORDER);

    printf("Reordering vectors...\n");
    spx_perm_t *p = spx_mat_get_perm(B);
    spx_vec_reorder(perm_x, p);
    spx_vec_reorder(perm_y, p);

    printf("Running the SpMV kernel...\n");
    spx_matvec_mult(0.1, B, perm_x, 0.41, perm_y);

    printf("Result: ");
    spx_vec_inv_reorder(perm_y, p);
    /* spx_vec_print(perm_y); */

    /* Cleanup */
    spx_input_destroy(input);
    spx_mat_destroy(A);
    spx_partition_destroy(parts);
    spx_vec_destroy(x);
    spx_vec_destroy(y);
    spx_vec_destroy(perm_x);
    spx_vec_destroy(perm_y);

    return 0;
}
