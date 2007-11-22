#include "matrix.h"
#include "sparse.h"
#include "method.h"
#include "spm_crsr.h"

static void spm_crs_multiply_ur(spm_t *matrix, vector_t *in, vector_t *out)
{
        elem_t *y = out->elements;
        elem_t *x = in->elements;
        elem_t *values = matrix->crs.values;
        index_t *row_ptr = matrix->crs.row_ptr;
        index_t *col_ind = matrix->crs.col_ind;
        unsigned long n = matrix->crs.nrows;
	unsigned long i0=0, i1=0, i2=0, i3=0, i4=0, i5=0, i6=0, i7=0;

        unsigned long i, j;

        __asm__ __volatile__(" # spmxv start");
        for(i=0; i<n; i++) {
		j = row_ptr[i];
                unsigned long j_end = row_ptr[i+1];
                unsigned long j_end_aligned = j + ((j_end - j) & (~(8-1)));

		//printf("j_end_aligned = %lu\n", j_end_aligned);
		__asm__ __volatile__(
			" movlpd (%[elem]), %%xmm8\n\t"
			" cmpq %[j_ea], %[j]\n\t"
			" jz .lead_out_ur \n\t"
			".beg_loop_ur: \n\t"
			" movq 0(%[col_ind]), %[i0]\n\t"
			" movq 8(%[col_ind]), %[i1]\n\t"
			" movq 16(%[col_ind]), %[i2]\n\t"
			" movq 24(%[col_ind]), %[i3]\n\t"
			" movq 32(%[col_ind]), %[i4]\n\t"
			" movq 40(%[col_ind]), %[i5]\n\t"
			" movq 48(%[col_ind]), %[i6]\n\t"
			" movq 56(%[col_ind]), %[i7]\n\t"
			" movlpd 0(%[vals]), %%xmm0\n\t"
			" movlpd 8(%[vals]), %%xmm1\n\t"
			" movlpd 16(%[vals]), %%xmm2\n\t"
			" movlpd 24(%[vals]), %%xmm3\n\t"
			" movlpd 32(%[vals]), %%xmm4\n\t"
			" movlpd 40(%[vals]), %%xmm5\n\t"
			" movlpd 48(%[vals]), %%xmm6\n\t"
			" movlpd 56(%[vals]), %%xmm7\n\t"
			" mulsd (%[x], %[i0], 8), %%xmm0\n\t"
			" mulsd (%[x], %[i1], 8), %%xmm1\n\t"
			" mulsd (%[x], %[i2], 8), %%xmm2\n\t"
			" mulsd (%[x], %[i3], 8), %%xmm3\n\t"
			" mulsd (%[x], %[i4], 8), %%xmm4\n\t"
			" mulsd (%[x], %[i5], 8), %%xmm5\n\t"
			" mulsd (%[x], %[i6], 8), %%xmm6\n\t"
			" mulsd (%[x], %[i7], 8), %%xmm7\n\t"
			" addsd %%xmm0, %%xmm8\n\t"
			" addsd %%xmm1, %%xmm8\n\t"
			" addsd %%xmm2, %%xmm8\n\t"
			" addsd %%xmm3, %%xmm8\n\t"
			" addsd %%xmm4, %%xmm8\n\t"
			" addsd %%xmm5, %%xmm8\n\t"
			" addsd %%xmm6, %%xmm8\n\t"
			" addsd %%xmm7, %%xmm8\n\t"
			" addq $8, %[j] \n\t"
			" addq $64, %[vals] \n\t"
			" addq $64, %[col_ind] \n\t"
			" cmpq %[j], %[j_ea]\n\t"
			" ja .beg_loop_ur\n\t"
			" .lead_out_ur:\n\t"
			" cmpq %[j], %[j_e]\n\t"
			" jz .end_ur\n\t"
			" movq 0(%[col_ind]), %[i0]\n\t"
			" movlpd 0(%[vals]), %%xmm0\n\t"
			" mulsd (%[x], %[i0], 8), %%xmm0\n\t"
			" addsd %%xmm0, %%xmm8\n\t"
			" addq $8, %[vals] \n\t"
			" addq $8, %[col_ind] \n\t"
			" incq  %[j] \n\t"
			" jmp .lead_out_ur\n\t"
			" .end_ur:\n\t"
			" movlpd %%xmm8, (%[elem])\n\t"
			: [j]  "+r" (j),
			  [i0] "=&r" (i0),
			  [i1] "=&r" (i1),
			  [i2] "=&r" (i2),
			  [i3] "=&r" (i3),
			  [i4] "=&r" (i4),
			  [i5] "=&r" (i5),
			  [i6] "=&r" (i6),
			  [i7] "=&r" (i7)
			: "[j]" (j),
			  [vals]    "r" (&values[j]),
			  [x]       "r" (x),
			  [j_ea]    "r" (j_end_aligned),
			  [j_e]     "r" (j_end),
			  [col_ind] "r" (&col_ind[j]),
			  [elem]    "r" (&y[i])
			: "xmm0", "xmm1", "xmm2", "xmm3",
			  "xmm4", "xmm5", "xmm6", "xmm7",
			  "xmm8",
			  "memory"
		);

        }
        __asm__ __volatile__(" # spmxv end");
}
METHOD_INIT(spm_crs_multiply_ur, spm_crs_load)

static void spm_crs_multiply_ur_vec(spm_t *matrix, vector_t *in, vector_t *out)
{
        elem_t *y = out->elements;
        elem_t *x = in->elements;
        elem_t *values = matrix->crs.values;
        index_t *row_ptr = matrix->crs.row_ptr;
        index_t *col_ind = matrix->crs.col_ind;
        unsigned long n = matrix->crs.nrows;
	unsigned long i0=0, i1=0, i2=0, i3=0, i4=0, i5=0, i6=0, i7=0, i8=0;

        unsigned long i, j;

        __asm__ __volatile__(" # spmxv start");
        for(i=0; i<n; i++) {
		j = row_ptr[i];
                unsigned long j_end = row_ptr[i+1];
                unsigned long j_end_aligned = j + ((j_end - j) & (~(8-1)));

		//printf("j_end_aligned = %lu\n", j_end_aligned);
		__asm__ __volatile__(
			" movlpd (%[elem]), %%xmm8\n\t"
			" cmpq %[j_ea], %[j]\n\t"
			" jz .lead_out_ur_vec \n\t"
			".beg_loop_ur_vec: \n\t"

			" movq 0(%[col_ind]), %[i0]\n\t"
			" movq 8(%[col_ind]), %[i1]\n\t"
			" movq 16(%[col_ind]), %[i2]\n\t"
			" movq 24(%[col_ind]), %[i3]\n\t"
			" movq 32(%[col_ind]), %[i4]\n\t"
			" movq 40(%[col_ind]), %[i5]\n\t"
			" movq 48(%[col_ind]), %[i6]\n\t"
			" movq 56(%[col_ind]), %[i7]\n\t"

			" movlpd (%[x], %[i0], 8), %%xmm0 \n\t"
			" movhpd (%[x], %[i1], 8), %%xmm0 \n\t"
			" movlpd (%[x], %[i2], 8), %%xmm1 \n\t"
			" movhpd (%[x], %[i3], 8), %%xmm1 \n\t"
			" movlpd (%[x], %[i4], 8), %%xmm2 \n\t"
			" movhpd (%[x], %[i5], 8), %%xmm2 \n\t"
			" movlpd (%[x], %[i6], 8), %%xmm3 \n\t"
			" movhpd (%[x], %[i7], 8), %%xmm3 \n\t"
			
			" mulpd 0(%[vals]), %%xmm0 \n\t"
			" mulpd 16(%[vals]), %%xmm1 \n\t"
			" mulpd 32(%[vals]), %%xmm2 \n\t"
			" mulpd 48(%[vals]), %%xmm3 \n\t"

			" haddpd %%xmm1, %%xmm0 \n\t"
			" haddpd %%xmm3, %%xmm2 \n\t"
			" haddpd %%xmm2, %%xmm0 \n\t"
			" haddpd %%xmm0, %%xmm0 \n\t" // first xmm0 is dummy
			" addsd  %%xmm0, %%xmm8 \n\t"

			" addq $8, %[j] \n\t"
			" addq $64, %[vals] \n\t"
			" addq $64, %[col_ind] \n\t"
			" cmpq %[j], %[j_ea]\n\t"
			" ja .beg_loop_ur_vec\n\t"
			" .lead_out_ur_vec:\n\t"
			" cmpq %[j], %[j_e]\n\t"
			" jz .end_ur_vec\n\t"
			" movq 0(%[col_ind]), %[i0]\n\t"
			" movlpd 0(%[vals]), %%xmm0\n\t"
			" mulsd (%[x], %[i0], 8), %%xmm0\n\t"
			" addsd %%xmm0, %%xmm8\n\t"
			" incq  %[j] \n\t"
			" addq $8, %[vals] \n\t"
			" addq $8, %[col_ind] \n\t"
			" jmp .lead_out_ur_vec\n\t"
			" .end_ur_vec:\n\t"
			" movlpd %%xmm8, (%[elem])\n\t"
			: [j]  "+r" (j),
			  [i0] "=&r" (i0),
			  [i1] "=&r" (i1),
			  [i2] "=&r" (i2),
			  [i3] "=&r" (i3),
			  [i4] "=&r" (i4),
			  [i5] "=&r" (i5),
			  [i6] "=&r" (i6),
			  [i7] "=&r" (i7)
			: "[j]" (j),
			  [vals]    "r" (&values[j]),
			  [x]       "r" (x),
			  [j_ea]    "r" (j_end_aligned),
			  [j_e]     "r" (j_end),
			  [col_ind] "r" (&col_ind[j]),
			  //[buffy]   "r" (buffy),
			  [elem]    "r" (&y[i])
			: "xmm0", "xmm1", "xmm2", "xmm3",
			  "xmm4", "xmm5", "xmm6", "xmm7",
			  "xmm8",
			  "memory"
		);

        }
        __asm__ __volatile__(" # spmxv end");
}
METHOD_INIT(spm_crs_multiply_ur_vec, spm_crs_load)

static void spm_crs_multiply_ur_vec2(spm_t *matrix, vector_t *in, vector_t *out)
{
        elem_t *y = out->elements;
        elem_t *x = in->elements;
        elem_t *values = matrix->crs.values;
        index_t *row_ptr = matrix->crs.row_ptr;
        index_t *col_ind = matrix->crs.col_ind;
        unsigned long n = matrix->crs.nrows;
	unsigned long i0=0, i1=0, i2=0, i3=0, i4=0, i5=0, i6=0, i7=0;

        unsigned long i, j;

        __asm__ __volatile__(" # spmxv start");
        for(i=0; i<n; i++) {
		j = row_ptr[i];
                unsigned long j_end = row_ptr[i+1];
                unsigned long j_end_aligned = j + ((j_end - j) & (~(8-1)));

		//printf("j_end_aligned = %lu\n", j_end_aligned);
		__asm__ __volatile__(
			" movlpd (%[elem]), %%xmm8\n\t"
			" cmpq %[j_ea], %[j]\n\t"
			" jz .lead_out_ur_vec2 \n\t"
			".beg_loop_ur_vec2: \n\t"

			" movq 0(%[col_ind]), %[i0]\n\t"
			" movlpd (%[x], %[i0], 8), %%xmm0 \n\t"
			" movq 8(%[col_ind]), %[i1]\n\t"
			" movhpd (%[x], %[i1], 8), %%xmm0 \n\t"
			" mulpd 0(%[vals]), %%xmm0 \n\t"
			" movq 16(%[col_ind]), %[i2]\n\t"
			" movlpd (%[x], %[i2], 8), %%xmm1 \n\t"
			" movq 24(%[col_ind]), %[i3]\n\t"
			" movhpd (%[x], %[i3], 8), %%xmm1 \n\t"
			" mulpd 16(%[vals]), %%xmm1 \n\t"
			" haddpd %%xmm1, %%xmm0 \n\t"
			" movq 32(%[col_ind]), %[i4]\n\t"
			" movlpd (%[x], %[i4], 8), %%xmm2 \n\t"
			" movq 40(%[col_ind]), %[i5]\n\t"
			" movhpd (%[x], %[i5], 8), %%xmm2 \n\t"
			" mulpd 32(%[vals]), %%xmm2 \n\t"
			" movq 48(%[col_ind]), %[i6]\n\t"
			" movlpd (%[x], %[i6], 8), %%xmm3 \n\t"
			" movq 56(%[col_ind]), %[i7]\n\t"
			" movhpd (%[x], %[i7], 8), %%xmm3 \n\t"
			" mulpd 48(%[vals]), %%xmm3 \n\t"
			" haddpd %%xmm3, %%xmm2 \n\t"
			" haddpd %%xmm2, %%xmm0 \n\t"
			" haddpd %%xmm0, %%xmm0 \n\t" // first xmm0 is dummy
			" addsd  %%xmm0, %%xmm8 \n\t"

			" addq $8, %[j] \n\t"
			" addq $64, %[vals] \n\t"
			" addq $64, %[col_ind] \n\t"
			" cmpq %[j], %[j_ea]\n\t"
			" ja .beg_loop_ur_vec2\n\t"
			" .lead_out_ur_vec2:\n\t"
			" cmpq %[j], %[j_e]\n\t"
			" jz .end_ur_vec2\n\t"
			" movq 0(%[col_ind]), %[i0]\n\t"
			" movlpd 0(%[vals]), %%xmm0\n\t"
			" mulsd (%[x], %[i0], 8), %%xmm0\n\t"
			" addsd %%xmm0, %%xmm8\n\t"
			" addq $8, %[vals] \n\t"
			" addq $8, %[col_ind] \n\t"
			" incq  %[j] \n\t"
			" jmp .lead_out_ur_vec2\n\t"
			" .end_ur_vec2:\n\t"
			" movlpd %%xmm8, (%[elem])\n\t"
			: [j]  "+r" (j),
			  [i0] "=&r" (i0),
			  [i1] "=&r" (i1),
			  [i2] "=&r" (i2),
			  [i3] "=&r" (i3),
			  [i4] "=&r" (i4),
			  [i5] "=&r" (i5),
			  [i6] "=&r" (i6),
			  [i7] "=&r" (i7)
			: "[j]" (j),
			  [vals]    "r" (&values[j]),
			  [x]       "r" (x),
			  [j_ea]    "r" (j_end_aligned),
			  [j_e]     "r" (j_end),
			  [col_ind] "r" (&col_ind[j]),
			  //[buffy]   "r" (buffy),
			  [elem]    "r" (&y[i])
			: "xmm0", "xmm1", "xmm2", "xmm3",
			  "xmm4", "xmm5", "xmm6", "xmm7",
			  "xmm8",
			  "memory"
		);

        }
        __asm__ __volatile__(" # spmxv end");
}
METHOD_INIT(spm_crs_multiply_ur_vec2, spm_crs_load)

static void spm_crs_multiply_ur_vec16(spm_t *matrix, vector_t *in, vector_t *out)
{
        elem_t *y = out->elements;
        elem_t *x = in->elements;
        elem_t *values = matrix->crs.values;
        index_t *row_ptr = matrix->crs.row_ptr;
        index_t *col_ind = matrix->crs.col_ind;
        unsigned long n = matrix->crs.nrows;
	unsigned long i0=0, i1=0, i2=0, i3=0, i4=0, i5=0, i6=0, i7=0, i8=0;

        unsigned long i, j;

        __asm__ __volatile__(" # spmxv start");
        for(i=0; i<n; i++) {
		j = row_ptr[i];
                unsigned long j_end = row_ptr[i+1];
                unsigned long j_end_aligned = j + ((j_end - j) & (~(16-1)));

		//printf("j=%lu j_end_aligned = %lu\n", j, j_end_aligned);
		__asm__ __volatile__(
			" movlpd (%[elem]), %%xmm8\n\t"
			" cmpq %[j_ea], %[j]\n\t"
			" jz .lead_out_ur_vec3 \n\t"
			".beg_loop_ur_vec3: \n\t"

			" movq 0(%[col_ind]), %[i0]\n\t"
			" movq 8(%[col_ind]), %[i1]\n\t"
			" movq 16(%[col_ind]), %[i2]\n\t"
			" movq 24(%[col_ind]), %[i3]\n\t"
			" movq 32(%[col_ind]), %[i4]\n\t"
			" movq 40(%[col_ind]), %[i5]\n\t"
			" movq 48(%[col_ind]), %[i6]\n\t"
			" movq 56(%[col_ind]), %[i7]\n\t"

			" movlpd (%[x], %[i0], 8), %%xmm0 \n\t"
			" movhpd (%[x], %[i1], 8), %%xmm0 \n\t"
			" movlpd (%[x], %[i2], 8), %%xmm1 \n\t"
			" movhpd (%[x], %[i3], 8), %%xmm1 \n\t"
			" movlpd (%[x], %[i4], 8), %%xmm2 \n\t"
			" movhpd (%[x], %[i5], 8), %%xmm2 \n\t"
			" movlpd (%[x], %[i6], 8), %%xmm3 \n\t"
			" movhpd (%[x], %[i7], 8), %%xmm3 \n\t"

			" movq 64(%[col_ind]), %[i0]\n\t"
			" movq 72(%[col_ind]), %[i1]\n\t"
			" movq 80(%[col_ind]), %[i2]\n\t"
			" movq 88(%[col_ind]), %[i3]\n\t"
			" movq 96(%[col_ind]), %[i4]\n\t"
			" movq 104(%[col_ind]), %[i5]\n\t"
			" movq 112(%[col_ind]), %[i6]\n\t"
			" movq 120(%[col_ind]), %[i7]\n\t"

			" movlpd (%[x], %[i0], 8), %%xmm4 \n\t"
			" movhpd (%[x], %[i1], 8), %%xmm4 \n\t"
			" movlpd (%[x], %[i2], 8), %%xmm5 \n\t"
			" movhpd (%[x], %[i3], 8), %%xmm5 \n\t"
			" movlpd (%[x], %[i4], 8), %%xmm6 \n\t"
			" movhpd (%[x], %[i5], 8), %%xmm6 \n\t"
			" movlpd (%[x], %[i6], 8), %%xmm7 \n\t"
			" movhpd (%[x], %[i7], 8), %%xmm7 \n\t"
			
			" mulpd 0(%[vals]), %%xmm0 \n\t"
			" mulpd 16(%[vals]), %%xmm1 \n\t"
			" mulpd 32(%[vals]), %%xmm2 \n\t"
			" mulpd 48(%[vals]), %%xmm3 \n\t"
			" mulpd 64(%[vals]), %%xmm4 \n\t"
			" mulpd 80(%[vals]), %%xmm5 \n\t"
			" mulpd 96(%[vals]), %%xmm6 \n\t"
			" mulpd 112(%[vals]), %%xmm7 \n\t"

			" haddpd %%xmm1, %%xmm0 \n\t"
			" haddpd %%xmm3, %%xmm2 \n\t"
			" haddpd %%xmm5, %%xmm4 \n\t"
			" haddpd %%xmm7, %%xmm6 \n\t"
			" haddpd %%xmm2, %%xmm0 \n\t"
			" haddpd %%xmm6, %%xmm4 \n\t"
			" haddpd %%xmm4, %%xmm0 \n\t"
			" haddpd %%xmm0, %%xmm0 \n\t" // first xmm0 is dummy
			" addsd  %%xmm0, %%xmm8 \n\t"

			" addq $16, %[j] \n\t"
			" addq $128, %[vals] \n\t"
			" addq $128, %[col_ind] \n\t"
			" cmpq %[j], %[j_ea]\n\t"
			" ja .beg_loop_ur_vec3\n\t"
			" .lead_out_ur_vec3:\n\t"
			" cmpq %[j], %[j_e]\n\t"
			" jz .end_ur_vec3\n\t"
			" movq 0(%[col_ind]), %[i0]\n\t"
			" movlpd 0(%[vals]), %%xmm0\n\t"
			" mulsd (%[x], %[i0], 8), %%xmm0\n\t"
			" addsd %%xmm0, %%xmm8\n\t"
			" addq $8, %[vals] \n\t"
			" addq $8, %[col_ind] \n\t"
			" incq  %[j] \n\t"
			" jmp .lead_out_ur_vec3\n\t"
			" .end_ur_vec3:\n\t"
			" movlpd %%xmm8, (%[elem])\n\t"
			: [j]  "+r" (j),
			  [i0] "=&r" (i0),
			  [i1] "=&r" (i1),
			  [i2] "=&r" (i2),
			  [i3] "=&r" (i3),
			  [i4] "=&r" (i4),
			  [i5] "=&r" (i5),
			  [i6] "=&r" (i6),
			  [i7] "=&r" (i7)
			: "[j]" (j),
			  [vals]    "r" (&values[j]),
			  [x]       "r" (x),
			  [j_ea]    "r" (j_end_aligned),
			  [j_e]     "r" (j_end),
			  [col_ind] "r" (&col_ind[j]),
			  //[buffy]   "r" (buffy),
			  [elem]    "r" (&y[i])
			: "xmm0", "xmm1", "xmm2", "xmm3",
			  "xmm4", "xmm5", "xmm6", "xmm7",
			  "xmm8",
			  "memory"
		);

        }
        __asm__ __volatile__(" # spmxv end");
}
METHOD_INIT(spm_crs_multiply_ur_vec16, spm_crs_load)


#if 0
void spm_crsr_multiply_urv(spm_crsr_t *crsr, vector_t *in, vector_t *out)
{
	elem_t *x = in->elements;
	elem_t *y = out->elements;
	const unsigned long nrows = crsr->rows;
	unsigned long r;

	for (r=0; r < nrows ; r++){
		const unsigned long units_nr = row->units_nr;
		spm_crsr_unit_t *unit = row->unit;
		yr = 0;
		const unsigned long unit_end = unit->len;
                const unsigned long unit_end_aligned = unit_end & (~(8-1));
		const elem_t 
	}
	row++;
	y[r] = yr;
}
#endif
