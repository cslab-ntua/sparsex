
#include <inttypes.h>

#include "dynarray.h"
#include "mmf.h"

#include "spm_csrdu.h"

typedef struct rle {
	uint64_t val;
	uint64_t freq;
} rle_t;

static void delta_encode(uint64_t *input, uint64_t *deltas, uint64_t size)
{
	uint64_t i;
	uint64_t prev = deltas[0] = input[0];
	for (i=1; i<size; i++){
		uint64_t curr = input[i];
		deltas[i] = curr - prev;
		prev = curr;
	}
}

static void rle_encode(uint64_t *in, uint64_t insize, rle_t *rles, uint64_t *rles_size)
{
	uint64_t rle_freq=1; // frequency of current value
	uint64_t rle_i=0;    // index of rle buffer
	rle_t *rle;
	uint64_t i;
	uint64_t prev = in[0];
	for (i=1; i<insize; i++){
		uint64_t curr = in[i];
		if (curr == prev){
			rle_freq++;
			continue;
		}

		rle = rles + rle_i++;
		rle->val = prev;
		rle->freq = rle_freq;
		rle_freq = 1;
		prev = curr;
	}

	rle = &rles[rle_i++];
	rle->val = in[insize-1];
	rle->freq = rle_freq;

	*rles_size = rle_i;
}

struct unit_state {
	uint64_t start, size; // unit start / size in the row
	uint8_t ci_size;      // current size (type) of unit
	char new_row;         // new row existance
	uint64_t *deltas;     // array of deltas for current row
};

static struct csrdu_st  {
	int sp_minlen; // min length for a SPARSE unit (0->use maximum size)
	int de_minlen; // min length for a DENSE unit (0->no dense units)
	struct unit_state unit;
	dynarray_t *da_ctl;
} state;

#define DE_MINLEN 0
#define SP_MINLEN 0
static void set_params()
{
	char *e;
	e = getenv("CSRDU_DE_MINLEN");
	state.de_minlen = e ? atoi(e) : DE_MINLEN;
	e = getenv("CSRDU_SP_MINLEN");
	state.sp_minlen = e ? atoi(e) : SP_MINLEN;

	printf("csrdu_params: sp_minlen:%d de_minlen:%d\n",state.sp_minlen,state.de_minlen);
}

static void de_add_unit()
{
	uint8_t *ctl_flags = dynarray_alloc_nr(state.da_ctl, 2);
	uint8_t *ctl_size = ctl_flags + 1;

	struct unit_state *ust = &state.unit;

	*ctl_flags = 0;
	*ctl_size = (uint8_t)ust->size;
	if (ust->new_row){
		 spm_csrdu_fl_setnr(ctl_flags);
		 ust->new_row = 0;
	}

	da_uc_put_ul(state.da_ctl, ust->deltas[ust->start]);
	ust->start += ust->size;
	ust->size = 0;
	ust->ci_size = SPM_CSRDU_CISIZE_U8;
}

static void sp_add_header()
{
	uint8_t *ctl_flags = dynarray_alloc_nr(state.da_ctl, 2);
	uint8_t *ctl_size = ctl_flags + 1;

	struct unit_state *ust = &state.unit;

	*ctl_flags = 0;
	spm_csrdu_fl_setsp(ctl_flags);
	spm_csrdu_fl_setcisize(ctl_flags, ust->ci_size);
	*ctl_size = (uint8_t)ust->size;
	if (ust->new_row){
		 spm_csrdu_fl_setnr(ctl_flags);
		 ust->new_row = 0;
	}
}

static void sp_add_body()
{
	struct unit_state *ust = &state.unit;
	uint64_t *src;
	void *dst;

	uint64_t ci_size = ust->ci_size;
	uint64_t dsize = spm_csrdu_cisize_bytes(ci_size);
	uint64_t usize = ust->size;

	dst = dynarray_alloc_nr(state.da_ctl, usize*dsize);
	src = ust->deltas + ust->start;

	spm_csrdu_cisize_copy(dst, src, usize, ci_size);
}

static void sp_add_unit()
{
	struct unit_state *ust = &state.unit;
	sp_add_header();
	sp_add_body();
	ust->start += ust->size;
	ust->size = 0;
	ust->ci_size = SPM_CSRDU_CISIZE_U8;
}


static void handle_row(uint64_t *deltas, uint64_t deltas_size,
                       rle_t *rles, uint64_t rles_size)
{
	uint64_t i=0;
	struct unit_state *ust = &state.unit;
	ust->start = 0;
	ust->size = 0;
	ust->ci_size = SPM_CSRDU_CISIZE_U8;
	ust->new_row = 1;
	ust->deltas = deltas;

	int sp_minlen = state.sp_minlen;
	int de_minlen = state.de_minlen;

	for (i=0; i < rles_size; i++){
		rle_t *rle = rles + i;

		// check if this is a large enough dense unit
		if (de_minlen && (rle->val == 1) && (rle->freq >= de_minlen)){
			// add previous sparse unit (if exists)
			if (ust->size){
				sp_add_unit();
			}
			// add dense unit
			ust->size = rle->freq;
			de_add_unit();
			continue;
		}

		// check if by adding this rle the usize gets > max unit size
		if (rle->freq + ust->size > SPM_CSRDU_SIZE_MAX){
			sp_add_unit();
		}

		// check if the ci_size changes with the new column index
		uint8_t new_ci_size = spm_csrdu_cisize(rle->val);
		if (new_ci_size > ust->ci_size) {
			// check the unit is large enough to be commited
			if (sp_minlen && (ust->size >= sp_minlen)){
				sp_add_unit();
			}
			ust->ci_size = new_ci_size;
		}

		ust->size += rle->freq;
	}

	// add remaining unit
	if (ust->size){
		sp_add_unit();
	}

}

SPM_CSRDU_TYPE *SPM_CSRDU_NAME(_init_mmf)(char *mmf_file,
                                          uint64_t *nrows, uint64_t *ncols, uint64_t *nnz)
{
	SPM_CSRDU_TYPE *csrdu;
	csrdu = malloc(sizeof(SPM_CSRDU_TYPE));
	if (!csrdu){
		perror("malloc");
		exit(1);
	}

	FILE *mmf = mmf_init(mmf_file, nrows, ncols, nnz);
	csrdu->nnz = *nnz;
	csrdu->ncols = *ncols;
	csrdu->nrows= *nrows;
	csrdu->values = malloc(sizeof(ELEM_TYPE)*csrdu->nnz);
	if (!csrdu->values){
		perror("malloc");
		exit(1);
	}

	dynarray_t *da_cis = dynarray_create(sizeof(uint64_t), 512);
	dynarray_t *da_deltas = dynarray_create(sizeof(uint64_t), 512);
	dynarray_t *da_rles = dynarray_create(sizeof(rle_t), 512);
	state.da_ctl = dynarray_create(sizeof(uint8_t), 4096);

	set_params();

	uint64_t row, col, row_prev, val_i=0;
	double val;
	row_prev = 0;
	//while (mmf_get_next(mmf, &row, &col, &val)) {
	for(;;){
		int ret = mmf_get_next(mmf, &row, &col, &val);
		if ((row != row_prev)|| !ret){
			uint64_t *cis = dynarray_get(da_cis, 0);
			uint64_t cis_size = dynarray_size(da_cis);
			uint64_t *deltas = dynarray_alloc_nr(da_deltas, cis_size);
			rle_t *rles = dynarray_alloc_nr(da_rles, cis_size);
			uint64_t rles_size;

			delta_encode(cis, deltas, cis_size);
			rle_encode(deltas, cis_size, rles, &rles_size);
			handle_row(deltas, cis_size, rles, rles_size);

			if (!ret){
				break;
			}

			dynarray_dealloc_all(da_cis);
			dynarray_dealloc_all(da_deltas);
			dynarray_dealloc_all(da_rles);
			row_prev = row;
		}

		uint64_t *ci = dynarray_alloc(da_cis);
		*ci = col;
		csrdu->values[val_i++] = (ELEM_TYPE)val;
	}

	free(dynarray_destroy(da_cis));
	free(dynarray_destroy(da_deltas));
	free(dynarray_destroy(da_rles));

	printf("ctl_size:%lu ", dynarray_size(state.da_ctl));
	csrdu->ctl = dynarray_destroy(state.da_ctl);
	printf("ctl: %p\n", csrdu->ctl);
	return csrdu;
}
