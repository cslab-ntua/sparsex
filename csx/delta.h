#ifndef CSX_DELTA_H__
#define CSX_DELTA_H__

namespace csx {

// Size for the delta-encoded units
typedef enum {
	DELTA_U8  = 0,
	DELTA_U16,
	DELTA_U32,
	DELTA_U64,
	DELTA_MAX
} DeltaSize;


#define  spm_csrdu_ucmax(ci_size) (1UL<<(8<<ci_size))
static inline DeltaSize getDeltaSize(uint64_t u)
{
	if ( u < spm_csrdu_ucmax(DELTA_U8)){
		return DELTA_U8;
	}

	if ( u < spm_csrdu_ucmax(DELTA_U16)){
		return DELTA_U16;
	}

	if ( u < spm_csrdu_ucmax(DELTA_U32)){
		return DELTA_U32;
	}

	return DELTA_U64;
}


// number of bytes required for specified DeltaSize
static inline int DeltaSize_getBytes(DeltaSize delta_size)
{
	switch (delta_size){
		case DELTA_U8:
		case DELTA_U16:
		case DELTA_U32:
		case DELTA_U64:
		return (1<<delta_size);
		default:
		assert(false);
	}
}

} // end csx namespace
#endif // CSX_DELTA_H__
