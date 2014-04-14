#include "sparsex.h"

int main(int argc, char **argv)
{
    double alpha = 0.8, beta = 0.42;
    int loops = 128;

    spx_init();

    /* Load matrix from MMF file */
    spx_input_t *input = spx_input_load_mmf(argv[1]);

    /* Transform to CSX with reordering enabled */
    spx_options_set_from_env();
    spx_matrix_t *A = spx_mat_tune(input, OP_REORDER);

    /* Create random x and y vectors */
    spx_partition_t *parts = spx_mat_get_parts(A);
    spx_vector_t *x = spx_vec_create_random(spx_mat_get_ncols(A), parts);
    spx_vector_t *y = spx_vec_create_random(spx_mat_get_nrows(A), parts);

    /* Reorder vectors */
    spx_perm_t *p = spx_mat_get_perm(A);
    spx_vec_reorder(x, p);
    spx_vec_reorder(y, p);

    /* Run 128 loops of the SpMV kernel */
    spx_timer_t t;
    double elapsed_time, flops;
    int i;
    
    spx_timer_clear(&t);
    spx_timer_start(&t);
    for (i = 0; i < loops; i++) {
        spx_matvec_mult(alpha, A, x, beta, y);
    }
    spx_timer_pause(&t);
    elapsed_time = spx_timer_get_secs(&t);
    flops = (double) (2*loops*spx_mat_get_nnz(A)) /
        ((double) 1000*1000*elapsed_time);

    printf("Elapsed time: %lf secs\n", elapsed_time);
    printf("FLOPS: %lf\n", flops);

    /* Restore original ordering of resulting vector */
    spx_vec_inv_reorder(y, p);

    /* Print result */
    /* spx_vec_print(y); */

    /* Cleanup */
    spx_input_destroy(input);
    spx_mat_destroy(A);
    spx_partition_destroy(parts);
    spx_vec_destroy(x);
    spx_vec_destroy(y);

    return 0;
}
