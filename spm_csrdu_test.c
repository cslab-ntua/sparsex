#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#include "spm_csrdu.h"

#define ALIGN(buf,a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))

int main(int argc, char **argv)
{
	if (argc < 2){
		fprintf(stderr, "Usage: %s <mmf_file>\n", argv[0]);
		exit(1);
	}

	uint64_t nrows0, ncols0, nnz0;
	spm_csrdu_double_t *csrdu;
	csrdu = spm_csrdu_double_init_mmf(argv[1], &nrows0, &ncols0, &nnz0);

	int align = 0;
	int jmp = 0;
	char *e;
	if ( (e = getenv("CSRDU_ALIGNED")) ){
		align = !!atoi(e);
	}
	if ( (e = getenv("CSRDU_JMP")) ){
		jmp = !!atoi(e);
	}
	//printf("params: align:%d jmp:%d\n", align, jmp);

	uint8_t *ctl = csrdu->ctl;
	uint64_t row=1, col=1;
	double *vals = csrdu->values;
	double *vals_end = vals + nnz0;
	uint64_t i=0;
	printf("%lu %lu %lu\n", nrows0, ncols0, nnz0);
	for (;;){
		//printf("new unit (%lu)\n", ctl - csrdu->ctl);
		uint8_t flags = u8_get(ctl);
		uint8_t size = u8_get(ctl);
		//printf("flags=%u size=%u\n", flags, size);

		if (spm_csrdu_fl_isnr(flags)){
			row++;
			col=1;
		}

		#define __fmt "%lu %lu %-14lf\n"
		#define ALIGN_CTL(align) do { ctl = ALIGN(ctl, align); } while(0)
		switch (flags & SPM_CSRDU_FL_UNIT_MASK){

			case SPM_CSRDU_FL_UNIT_DENSE:
			col += uc_get_ul(ctl);
			for (i=0; i<size; i++){
				printf(__fmt, row, col+i, *vals++);
			}
			col += (size -1);
			break;

			#define PRINT_CSRDU_SP(bits)                       \
			if (!jmp) {                                        \
				if (align)                                 \
					ALIGN_CTL(bits/8);                 \
				for (i=0; i<size; i++){                    \
					col += u ## bits ## _get(ctl);     \
					printf(__fmt, row, col, *vals++);  \
				}                                          \
			} else {                                           \
				col += uc_get_ul(ctl);                     \
				if (align)                                 \
					ALIGN_CTL(bits/8);                 \
				printf(__fmt, row, col, *vals++);          \
				for (i=1; i<size; i++){                    \
					col += u ## bits ## _get(ctl);     \
					printf(__fmt, row, col, *vals++);  \
				}                                          \
			}                                                  \

			case SPM_CSRDU_FL_UNIT_SP_U8:
			PRINT_CSRDU_SP(8)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U16:
			PRINT_CSRDU_SP(16)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U32:
			PRINT_CSRDU_SP(32)
			break;

			case SPM_CSRDU_FL_UNIT_SP_U64:
			PRINT_CSRDU_SP(64)
			break;

			#undef PRINT_CSRDU_SP
			default:
			assert(0);
		}
		#undef ALIGN_CTL


		if (vals == vals_end)
			break;
	}
	spm_csrdu_double_destroy(csrdu);

	return 0;
}
