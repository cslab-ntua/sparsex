#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#include "dynarray.h"
#include "mmf.h"

#include "spm_csrdu.h"
#define MIN(x,y) (x < y ? x : y)

typedef struct rle {
	uint64_t val;
	uint64_t freq;
} rle_t;

static void delta_encode(uint64_t *input, uint64_t *deltas, uint64_t size)
{
	uint64_t i;
	uint64_t prev = input[0];
	deltas[0] = prev ;//+ 1; // make the 0,1,2,3 case => 1,1,1,1
	for (i=1; i<size; i++){
		uint64_t curr = input[i];
		deltas[i] = curr - prev;
		prev = curr;
	}
	//printf("prev=%lu input=%lu\n", deltas[0], input[0]);
}

static void rle_encode(uint64_t *in, uint64_t insize, rle_t *rles, uint64_t *rles_size)
{
	uint64_t rle_freq=1; // frequency of current value
	uint64_t rle_i=0;    // index of rle buffer
	rle_t *rle;
	uint64_t i;
	uint64_t prev = in[0];
	//printf("prev=%lu\n", prev);
	for (i=1; i<insize; i++){
		uint64_t curr = in[i];
		if (curr == prev){
			rle_freq++;
			continue;
		}

		rle = rles + rle_i++;
		rle->val = prev;
		rle->freq = rle_freq;
		//printf("rles[%lu] = (f:%lu v:%lu)\n", rle-rles, rle->freq, rle->val);
		rle_freq = 1;
		prev = curr;
	}

	rle = &rles[rle_i++];
	rle->val = in[insize-1];
	rle->freq = rle_freq;
	//printf("rles[%lu] = (f:%lu v:%lu)\n", rle-rles, rle->freq, rle->val);

	*rles_size = rle_i;
}

struct unit_state {
	uint64_t start, size; // unit start / size in the row
	uint64_t jmp;
	uint8_t ci_size;      // current size (type) of unit
	char new_row;         // new row existance
	uint64_t *deltas;     // array of deltas for current row
};

static struct stats {
	uint64_t units_de;
	uint64_t units_sp[SPM_CSRDU_CISIZE_NR];
} stats = {0};

static struct csrdu_st  {
	int sp_minlen;      // min length for a SPARSE unit (0->use maximum size)
	int de_minlen;      // min length for a DENSE unit (0->no dense units)
	uint64_t row_size;  // size of current row
	struct unit_state unit;
	dynarray_t *da_ctl;

	unsigned aligned:1;    // use aligned deltas
	unsigned jmp:1;        // use jumps
	unsigned verbose:1;
} state;

// verbose message
#define vmsg(fmt,args...) do { if (state.verbose){ printf(fmt, args); }} while(0)

// parameter defaults
#define DE_MINLEN_DEF 0
#define SP_MINLEN_DEF 0
#define ALIGNED_DEF 0
#define JMP_DEF 0
#define VERBOSE_DEF 1
static void set_params()
{
	char *e;
	e = getenv("CSRDU_DE_MINLEN");
	state.de_minlen = e ? atoi(e) : DE_MINLEN_DEF;
	e = getenv("CSRDU_SP_MINLEN");
	state.sp_minlen = e ? atoi(e) : SP_MINLEN_DEF;
	e = getenv("CSRDU_ALIGNED");
	state.aligned = e ? !!atoi(e) : ALIGNED_DEF;
	e = getenv("CSRDU_JMP");
	state.jmp = e ? !!atoi(e) : JMP_DEF;
	e = getenv("CSRDU_VERBOSE");
	state.verbose = e ? !!atoi(e) : VERBOSE_DEF;

	vmsg("csrdu_params: sp_minlen:%d de_minlen:%d aligned:%d jmp:%d\n",
              state.sp_minlen, state.de_minlen, state.aligned, state.jmp);
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

	//printf("de_add usize: %lu ctl offset: %lu\n", ust->size, dynarray_size(state.da_ctl));
	da_uc_put_ul(state.da_ctl, ust->deltas[ust->start]);
	ust->start += ust->size;
	ust->size = 0;
	ust->ci_size = SPM_CSRDU_CISIZE_U8;

	stats.units_de++;
}

static void sp_add_header(uint64_t usize, uint8_t ci_size, char *new_row_ptr)
{
	//printf("sp_add_header usize: %lu ctl offset: %lu\n", usize, dynarray_size(state.da_ctl));
	uint8_t *ctl_flags = dynarray_alloc_nr(state.da_ctl, 2);
	uint8_t *ctl_size = ctl_flags + 1;

	*ctl_flags = 0;
	spm_csrdu_fl_setsp(ctl_flags);
	spm_csrdu_fl_setcisize(ctl_flags, ci_size);
	*ctl_size = (uint8_t)usize;
	if (*new_row_ptr){
		 spm_csrdu_fl_setnr(ctl_flags);
		 *new_row_ptr = 0;
	}

	stats.units_sp[ci_size]++;
}

static void sp_add_body(uint64_t ustart, uint64_t usize, uint8_t ci_size, uint64_t *deltas)
{
	uint64_t *src;
	void *dst;

	uint64_t dsize = spm_csrdu_cisize_bytes(ci_size);

	dst = (state.aligned) ?
	      dynarray_alloc_nr_aligned(state.da_ctl, usize*dsize,  dsize):
	      dynarray_alloc_nr(state.da_ctl, usize*dsize);
	src = deltas + ustart;

	/*
	printf("sp_add_body: adding deltas: ");
	int i;
	for (i=0; i<usize; i++){
		printf("%lu ", src[i]);
	}
	printf("\n");
	*/

	spm_csrdu_cisize_copy(dst, src, usize, ci_size);
}

static void sp_add_unit()
{
	struct unit_state *ust = &state.unit;
	sp_add_header(ust->size, ust->ci_size, &ust->new_row);
	sp_add_body(ust->start, ust->size, ust->ci_size, ust->deltas);
	ust->start += ust->size;
	ust->size = 0;
	ust->ci_size = SPM_CSRDU_CISIZE_U8;
}

static void sp_jmp_add_unit()
{
	struct unit_state *ust = &state.unit;
	//printf("sp_jmp_add: ust->start: %lu ust->size: %lu state.row_size:%lu new_r:%d ctl_off:%lu\n", ust->start, ust->size, state.row_size, ust->new_row, dynarray_size(state.da_ctl));
	assert(ust->start + ust->size <= state.row_size);

	assert(ust->size > 0);
	uint64_t ustart = ust->start;
	uint64_t usize = ust->size;
	uint64_t ci_size = ust->ci_size;
	sp_add_header(usize, ci_size, &ust->new_row);
	da_uc_put_ul(state.da_ctl, ust->jmp);
	if (usize > 1)
		sp_add_body(ustart +1, usize -1, ci_size, ust->deltas);
	ust->start += ust->size;
	ust->size = 0;
	ust->ci_size = SPM_CSRDU_CISIZE_U8;
}
typedef void void_fn_t(void);
static void handle_row(uint64_t *deltas, uint64_t deltas_size,
                       rle_t *rles, uint64_t rles_size)
{
	struct unit_state *ust = &state.unit;
	ust->start = 0;
	ust->size = 0;
	ust->ci_size = SPM_CSRDU_CISIZE_U8;
	ust->new_row = 1;
	ust->deltas = deltas;

	int sp_minlen = state.sp_minlen;
	int de_minlen = state.de_minlen;

	void_fn_t *sp_add = state.jmp ? sp_jmp_add_unit : sp_add_unit;

	// when state.jmp is true, then we have the following semantics
	rle_t *rle = rles;
	rle_t *rle_end = rles + rles_size;
	if (state.jmp){
		// first element is a jmp
		ust->jmp = ust->deltas[0];
		assert(rle->freq >= 1);
		if (rle->freq == 1){
			rle++;
		} else {
			rle->freq--;
		}
		ust->size = 1;
	}

	while (rle < rle_end) {
		//printf("rle->freq %lu rle->val %lu ust->size %lu\n", rle->freq, rle->val, ust->size);
		// check if this is a large enough dense unit
		if (de_minlen && (rle->val == 1) && (rle->freq >= (de_minlen-1))){
			// add previous sparse unit (if exists)
			if (ust->size > 1){
				// the next unit is dense, so there is always a jmp
				ust->size--;
				sp_add();
				ust->size = 1;
			}

			// add dense unit
			do {
				uint64_t s = MIN(rle->freq, SPM_CSRDU_SIZE_MAX - ust->size);
				ust->size += s;
				rle->freq -= s;
				de_add_unit();
			} while (rle->freq >= de_minlen);

			if (rle->freq == 0){
				rle++;
				if (rle == rle_end){
					break;
				}
			}

			if (state.jmp){
				rle->freq--;
				ust->jmp = ust->deltas[ust->start];
				ust->size = 1;
				if (rle == rle_end){
					break;
				}
				continue;
			}
		}

		// check if the ci_size changes with the new column index
		uint8_t new_ci_size = spm_csrdu_cisize(rle->val);
		if (new_ci_size > ust->ci_size) {
			// check the unit is large enough to be commited
			if (sp_minlen && (ust->size >= sp_minlen)){
				sp_add();
				assert(0);
			}
			ust->ci_size = new_ci_size;
		}


		// check if by adding this rle the usize gets > max unit size
		while (rle->freq + ust->size > SPM_CSRDU_SIZE_MAX) {
			rle->freq -= (SPM_CSRDU_SIZE_MAX - ust->size);
			ust->size = SPM_CSRDU_SIZE_MAX;
			sp_add();
			if (state.jmp){
				ust->jmp = ust->deltas[ust->start];
				ust->size = 1;
				rle->freq--;
			}
		}

		ust->size += rle->freq;
		rle++;
	}

	// add remaining unit
	if (ust->size){
		sp_add();
	}
}

void SPM_CSRDU_NAME(_destroy)(SPM_CSRDU_TYPE *csrdu)
{
	free(csrdu->values);
	free(csrdu->ctl);
	free(csrdu);
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
	for(;;){
		int ret = mmf_get_next(mmf, &row, &col, &val);
		if ((row != row_prev)|| !ret){
			uint64_t *cis = dynarray_get(da_cis, 0);
			uint64_t cis_size = dynarray_size(da_cis);
			uint64_t *deltas = dynarray_alloc_nr(da_deltas, cis_size);
			rle_t *rles = dynarray_alloc_nr(da_rles, cis_size);
			uint64_t rles_size;

			state.row_size = cis_size;
			delta_encode(cis, deltas, cis_size);
			rle_encode(deltas, cis_size, rles, &rles_size);
			//printf("--- handle row start (row=%lu,row_size=%lu)\n", row_prev, state.row_size);
			handle_row(deltas, cis_size, rles, rles_size);
			//printf("--- handle row end\n");

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

	fclose(mmf);
	free(dynarray_destroy(da_cis));
	free(dynarray_destroy(da_deltas));
	free(dynarray_destroy(da_rles));
	assert(val_i == csrdu->nnz);

	vmsg("ctl_size: %lu \n", dynarray_size(state.da_ctl));
	vmsg("units:\n de    :%6lu\n sp(8) :%6lu\n sp(16):%6lu\n sp(32):%6lu\n sp(64):%6lu\n",
	      stats.units_de, stats.units_sp[SPM_CSRDU_CISIZE_U8], stats.units_sp[SPM_CSRDU_CISIZE_U16], stats.units_sp[SPM_CSRDU_CISIZE_U32], stats.units_sp[SPM_CSRDU_CISIZE_U64]);
	csrdu->ctl = dynarray_destroy(state.da_ctl);
	return csrdu;
}
