#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

#include "phash.h"
#include "dynarray.h"
#include "spm_crs_vi.h"

void SPM_VI_NAME(_vals_to_vis)(ELEM_TYPE *vals, uint64_t vals_size,
                               SPM_CRSVI_VI_TYPE **uvals_idxs, ELEM_TYPE **uvals,
                               uint64_t *uvals_size)
{
	phash_t *vhash; /* hash for unique values */
	dynarray_t *da_uvals; /* dynamic array for unique values */
	SPM_CRSVI_VI_TYPE *idxs; /* value indices */
	unsigned long key, idx;
	uint64_t uvals_nr, i;

	da_uvals = dynarray_create(sizeof(ELEM_TYPE), 4096);
	assert(sizeof(ELEM_TYPE) <= sizeof(unsigned long));
	vhash = phash_new(12);
	idxs = malloc(sizeof(SPM_CRSVI_VI_TYPE)*vals_size);
	if (!idxs){
		perror("malloc");
		exit(1);
	}

	for (uvals_nr=0,i=0; i<vals_size; i++){
		key =  0;
		ELEM_TYPE val = vals[i];
		memcpy(&key, &val, sizeof(val));
		if (!phash_lookup(vhash, key, &idx)){ // need to insert value
			idx = phash_elements(vhash);
			phash_insert(vhash, key, idx);
			assert(idx == dynarray_size(da_uvals));
			assert(idx <= (1UL<<SPM_CRSVI_VI_BITS));
			ELEM_TYPE *uval = dynarray_alloc(da_uvals);
			*uval = val;
			uvals_nr++;
		 }
		 idxs[i] = (SPM_CRSVI_VI_TYPE)idx;
	}

	*uvals_size = dynarray_size(da_uvals);
	*uvals = dynarray_destroy(da_uvals);
	*uvals_idxs = idxs;
	phash_free(vhash);
}
