#include <map>
#include <algorithm>

#include <boost/foreach.hpp>

#define FOREACH BOOST_FOREACH

extern "C" {
	#include "dynarray.h"
}

#include "spm.h"
#include "delta.h"
#include "ctl.h"

using namespace csx;

template <typename IterT, typename ValT>
void DeltaEncode(IterT start, IterT end, ValT &x0)
{
	IterT i;
	ValT prev, tmp;

	prev = x0;
	for (i=start; i != end; ++i){
		tmp = *i;
		*i -= prev;
		prev = tmp;
	}
}

template <typename T>
void Copy(T *dst, uint64_t *src, long nr_items)
{
	for (long i=0; i<nr_items; i++){
		dst[i] = static_cast<T>(src[i]);
	}
}


#define LONGUC_SHIFT (7)
static inline void da_put_ul(dynarray_t *da, unsigned long val)
{
	uint8_t *uc;
	const unsigned shift = LONGUC_SHIFT;

	for (;;) {
		uc = (uint8_t *)dynarray_alloc(da);
		*uc = (val & ((1<<shift) - 1));
		if ( val < (1<<shift) ){
			break;
		}
		*uc |= (1<<shift);
		val >>= shift;
	}
}

#define u8_get(ptr) ({    \
	uint8_t *_ptr = ptr;  \
	ptr++; *_ptr;         \
})

#define u16_get(ptr) ({                 \
	uint16_t ret = *((uint16_t *)ptr);  \
	ptr += sizeof(uint16_t);            \
	ret;                                \
})

#define u32_get(ptr) ({                 \
	uint32_t ret = *((uint32_t *)ptr);  \
	ptr += sizeof(uint32_t);            \
	ret;                                \
})

#define u64_get(ptr) ({                 \
	uint64_t ret = *((uint64_t *)ptr);  \
	ptr += sizeof(uint64_t);            \
	ret;                                \
})

#define uc_get_ul(ptr)                  \
({                                      \
	unsigned long _val;                 \
	                                    \
	_val = u8_get(ptr);                 \
	if ( _val > 127 ){                  \
		unsigned shift = 7;             \
		unsigned long _uc;              \
		_val -= 128;                    \
		for (;;){                       \
			_uc = u8_get(ptr);          \
			if ( _uc > 127 ){           \
				_uc -= 128;             \
				_val += (_uc<<shift);   \
				shift += 7;             \
			} else {                    \
				_val += (_uc<<shift);   \
				break;                  \
			}                           \
		}                               \
	}                                   \
	_val;                               \
})                                      \

uint8_t CtlManager::getFlag(long pattern_id, uint64_t nnz)
{
	CtlManager::PatMap::iterator pi;
	uint8_t ret;
	pi = this->patterns.find(pattern_id);
	if (pi == this->patterns.end()){
		ret = this->flag_avail++;
		CtlManager::PatInfo patinfo(ret, nnz);
		this->patterns[pattern_id] = patinfo;
	} else {
		ret = pi->second.flag;
		pi->second.nr += nnz;
	}

	return ret;
}

//******** Static Pattern Id Mapping:
//  8  -> delta 8
//  16 -> delta 16
//  32 -> delta 32
//  64 -> delta 64
#define PID_DELTA_BASE 0

uint8_t *CtlManager::mkCtl()
{
	uint8_t *ret;
	this->ctl_da = dynarray_create(sizeof(uint8_t), 512);
	FOREACH(SpmRowElems &row, this->spm->rows){
		if (row.empty()){
			this->empty_rows++;
			continue;
		}
		this->doRow(row);
		row.clear();
	}

	ret = (uint8_t *)dynarray_destroy(this->ctl_da);
	return ret;
}

// Note that this function may allocate space in ctl_da
void CtlManager::updateNewRow(uint8_t *flags)
{
	if (!this->new_row)
		return;

	set_bit(flags, CTL_NR_BIT);
	this->new_row = false;
	if (this->empty_rows != 0){
		set_bit(flags, CTL_RJMP_BIT);
		da_put_ul(this->ctl_da, this->empty_rows);
		this->empty_rows = 0;
	}
}

void CtlManager::AddXs(std::vector<uint64_t> xs)
{
	uint8_t *ctl_flags, *ctl_size;
	long pat_id, xs_size, delta_bytes;
	uint64_t last_col, max;
	DeltaSize delta_size;
	std::vector<uint64_t>::iterator vi;
	void *dst;

	// do delta encoding
	xs_size = xs.size();
	last_col = xs[xs_size - 1];
	DeltaEncode(xs.begin(), xs.end(), this->last_col);
	this->last_col = last_col;

	// calculate the delta's size and the pattern id
	max = 0;
	if (xs_size > 1){
		vi = xs.begin();
		std::advance(vi, 1); // advance over jmp
		max = *(std::max_element(vi, xs.end()));
	}
	delta_size =  getDeltaSize(max);
	pat_id = delta_size + PID_DELTA_BASE;

	// set flags
	ctl_flags = (uint8_t *)dynarray_alloc_nr(this->ctl_da, 2);
	*ctl_flags = this->getFlag(PID_DELTA_BASE + pat_id, xs_size);
	this->updateNewRow(ctl_flags);

	// set size
	ctl_size = ctl_flags + 1;
	assert( (xs_size > 0) && (xs_size <= CTL_SIZE_MAX));
	*ctl_size = xs_size;

	// add jmp and deltas
	da_put_ul(this->ctl_da, xs[0]);

	//add deltas (if needed)
	if (xs_size > 1){
		delta_bytes = DeltaSize_getBytes(delta_size);
		dst = dynarray_alloc_nr_aligned(this->ctl_da, delta_bytes*(xs_size-1), delta_bytes);
		switch (delta_size){
			case DELTA_U8:  Copy((uint8_t  *)dst, &xs[1], xs_size-1); break;
			case DELTA_U16: Copy((uint16_t *)dst, &xs[1], xs_size-1); break;
			case DELTA_U32: Copy((uint32_t *)dst, &xs[1], xs_size-1); break;
			default:assert(false);
		}
	}

	xs.clear();
	return;
}

void CtlManager::AddPattern(const SpmRowElem &elem, uint64_t jmp)
{
	uint8_t *ctl_flags, *ctl_size;
	long pat_size, pat_id;
	uint64_t ujmp;

	pat_size = elem.pattern->getSize();
	pat_id = elem.pattern->getPatId();

	ctl_flags = (uint8_t *)dynarray_alloc_nr(this->ctl_da, 2);
	*ctl_flags = this->getFlag(pat_id, pat_size);
	this->updateNewRow(ctl_flags);

	ctl_size = ctl_flags + 1;
	assert(pat_size + 1 <= CTL_SIZE_MAX);
	// if there is a jmp we are implicitly including one more element
	// see also: 1feee866421a129fae861f094f64b6d803ecb8d5
	*ctl_size = pat_size + (jmp ? 1 : 0);

	assert(elem.x > this->last_col);
	ujmp = elem.x - this->last_col;
	da_put_ul(this->ctl_da, ujmp);
}

// Ctl Rules
// 1. Each unit leaves the x index at the last element it calculated on the
// current row
// 2. Size is the number of elements taht will be calculated
void CtlManager::doRow(const SpmRowElems &row)
{
	std::vector<uint64_t> xs;
	this->new_row = true;
	this->last_col = 0;
	FOREACH(const SpmRowElem &spm_elem, row){
		// check if this element contains a pattern
		if (spm_elem.pattern != NULL){
			uint64_t jmp;
			if (xs.size() > 0){
				jmp = xs.back();
				xs.pop_back();
			} else {
				jmp = 0;
			}
			if (xs.size() > 0)
				this->AddXs(xs);
			this->AddPattern(spm_elem, jmp);
			continue;
		}

		// check if we exceeded the maximum size for a unit
		assert(xs.size() <= CTL_SIZE_MAX);
		if (xs.size() == CTL_SIZE_MAX){
			this->AddXs(xs);
			continue;
		}

		xs.push_back(spm_elem.x);
	}

	if (xs.size() > 0)
		this->AddXs(xs);
}
