#include <sparsex/sparsex.h>

int main(int argc, char **argv)
{
    spx_value_t alpha = 0.8, beta = 0.42;
    const size_t nr_loops = 128;

    spx_init();

    /* Load matrix from MMF file */
    spx_input_t *input = spx_input_load_mmf(argv[1]);

    /* Transform to CSX */
    spx_option_set("spx.rt.nr_threads", "1");
    spx_option_set("spx.rt.cpu_affinity", "0");
    spx_option_set("spx.preproc.xform", "all");
    spx_option_set("spx.preproc.sampling", "portion");
    spx_option_set("spx.preproc.sampling.nr_samples", "48");
    spx_option_set("spx.preproc.sampling.portion", "0.01");
    /* spx_option_set("spx.matrix.symmetric", "true"); */
    spx_matrix_t *A = spx_mat_tune(input);

    /* Create random x and y vectors */
    spx_partition_t *parts = spx_mat_get_partition(A);
    spx_vector_t *x = spx_vec_create_random(spx_mat_get_ncols(A), parts);
    spx_vector_t *y = spx_vec_create_random(spx_mat_get_nrows(A), parts);

    /* Run 128 loops of the SpMV kernel */
    spx_timer_t t;
    double elapsed_time, flops;
    size_t i;

    spx_timer_clear(&t);
    spx_timer_start(&t);
    for (i = 0; i < nr_loops; i++) {
        spx_matvec_kernel(alpha, A, x, beta, y);
    }
    spx_timer_pause(&t);

    elapsed_time = spx_timer_get_secs(&t);
    flops = (double) (2 * nr_loops * spx_mat_get_nnz(A)) /
        ((double) 1000 * 1000 * elapsed_time);

    printf("Elapsed time: %lf secs\n", elapsed_time);
    printf("FLOPS: %lf\n", flops);

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
