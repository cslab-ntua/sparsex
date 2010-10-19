#include <map>
#include <algorithm>

#include <boost/foreach.hpp>

#define FOREACH BOOST_FOREACH

extern "C" {
	#include "dynarray.h"
}

#include "spm.h"
#include "delta.h"
#include "csx.h"

using namespace csx;

static bool debug = false;

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

uint8_t CsxManager::getFlag(long pattern_id, uint64_t nnz)
{
    CsxManager::PatMap::iterator pi;
    uint8_t ret;

    pi = this->patterns.find(pattern_id);
    if (pi == this->patterns.end()){
        ret = this->flag_avail++;
        CsxManager::PatInfo patinfo(ret, nnz);
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

csx_double_t *CsxManager::mkCsx()
{
    csx_double_t *csx;
    SPM *Spm;

    Spm = this->spm;
    csx = (csx_double_t *)malloc(sizeof(csx_double_t));
    this->values = (double *)malloc(sizeof(double)*Spm->nnz);
    if (!csx || !this->values){
        perror("malloc");
        exit(1);
    }
    this->ctl_da = dynarray_create(sizeof(uint8_t), 512);

    csx->nnz = Spm->nnz;
    csx->nrows = Spm->nrows;
    csx->ncols = Spm->ncols;
    csx->row_start = Spm->row_start;

    this->values_idx = 0;
    this->new_row = false;		// do not mark first row
    for (uint64_t i=0; i < Spm->getNrRows(); i++){
        const SpmRowElem *rbegin, *rend;

        rbegin = Spm->rbegin(i);
        rend = Spm->rend(i);
        if (debug)
            std::cerr << "mkCsx(): row: " << i << "\n";
        if (rbegin == rend){ 		// check if row is empty
            if (debug)
                std::cerr << "mkCsx(): row is empty" << std::endl;
            if (this->new_row == false){
                this->new_row = true; 	// in case the first row is empty
            } else {
                this->empty_rows++;
            }
            continue;
        }
        this->doRow(rbegin, rend);
        this->new_row = true;
    }
    csx->ctl_size = dynarray_size(this->ctl_da);
    //std::cerr << "csx->ctl_size=" << csx->ctl_size << "\n";
    csx->ctl = (uint8_t *)dynarray_destroy(this->ctl_da);
    //std::cerr << "csx->ctl=" << (unsigned long)csx->ctl << "\n";
    this->ctl_da = NULL;
    assert(this->values_idx == Spm->nnz);
    csx->values = this->values;
    this->values = NULL;
    this->values_idx = 0;

    return csx;
}

// Note that this function may allocate space in ctl_da
void CsxManager::updateNewRow(uint8_t *flags)
{
	if (!this->new_row)
		return;

	set_bit(flags, CTL_NR_BIT);
	this->new_row = false;
	if (this->empty_rows != 0){
		set_bit(flags, CTL_RJMP_BIT);
		da_put_ul(this->ctl_da, this->empty_rows + 1);
		this->empty_rows = 0;
		this->row_jmps = true;
	}
}

void CsxManager::AddXs(std::vector<uint64_t> &xs)
{
    //std::cerr << "AddXs: size:" << xs.size() << " first:" << xs[0] << "\n";
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
        std::advance(vi, 1); 	// advance over jmp
        max = *(std::max_element(vi, xs.end()));
    }
    delta_size =  getDeltaSize(max);
    pat_id = (8<<delta_size) + PID_DELTA_BASE;

    // set flags
    ctl_flags = (uint8_t *)dynarray_alloc_nr(this->ctl_da, 2);
    *ctl_flags = this->getFlag(PID_DELTA_BASE + pat_id, xs_size);

    // set size
    ctl_size = ctl_flags + 1;
    assert( (xs_size > 0) && (xs_size <= CTL_SIZE_MAX));
    *ctl_size = xs_size;

    // ctls_size, ctl_flags are not valid after this call
    this->updateNewRow(ctl_flags);

    // add jmp and deltas
    //std::cerr << "AddXs: Delta jmp:" << xs[0] << "\n";
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

void CsxManager::AddPattern(const SpmRowElem &elem, uint64_t jmp)
{
    uint8_t *ctl_flags, *ctl_size;
    long pat_size, pat_id;
    uint64_t ujmp;

    pat_size = elem.pattern->getSize();
    if (debug)
        std::cerr << "AddPattern jmp: " << jmp << " pat_size: " << pat_size << "\n";
    pat_id = elem.pattern->getPatId();
    ctl_flags = (uint8_t *)dynarray_alloc_nr(this->ctl_da, 2);
    *ctl_flags = this->getFlag(pat_id, pat_size);
    ctl_size = ctl_flags + 1;
    // std::cerr << pat_size << std::endl;
    assert(pat_size + (jmp ? 1 : 0) <= CTL_SIZE_MAX);
    // if there is a jmp we are implicitly including one more element
    // see also: 1feee866421a129fae861f094f64b6d803ecb8d5
    *ctl_size = pat_size + (jmp ? 1 : 0);
    // ctl_flags and ctl_size are not valid after this call
    this->updateNewRow(ctl_flags);
    //ujmp = elem.x - this->last_col;
    ujmp = jmp ? jmp : elem.x - this->last_col;
    if (debug)
        std::cerr << "AddPattern ujmp " << ujmp << "\n";
    da_put_ul(this->ctl_da, ujmp);
    this->last_col = elem.pattern->x_increase_jmp(this->spm->type, elem.x);
    if (debug)
        std::cerr << "last_col:" << this->last_col << "\n";
}

// return ujmp
uint64_t CsxManager::PreparePat(std::vector<uint64_t> &xs, const SpmRowElem &elem)
{
    uint64_t lastx;
    if (xs.size() == 0)
        return 0;

    if (elem.pattern->type != this->spm->type){
        this->AddXs(xs);
        return 0;
    }
    lastx = xs.back();
    // normaly we wouldn't need to check for this, since
    // it is assured by the parsing. Nevertheless, the
    // previous element can ``disappear'' if it is included
    // in another type of pattern.
    // Todo: maybe it's cleaner to fix the parsing
    if (elem.pattern->getNextX(lastx) != elem.x){
        this->AddXs(xs);
        return 0;
    }
    xs.pop_back();
    if (xs.size() > 0)
        this->AddXs(xs);
    return lastx - this->last_col;
}

// Ctl Rules
// 1. Each unit leaves the x index at the last element it calculated on the
// current row
// 2. Size is the number of elements that will be calculated
void CsxManager::doRow(const SpmRowElem *rbegin, const SpmRowElem *rend)
{
    std::vector<uint64_t> xs;
    uint64_t jmp;

    this->last_col = 1;
    for (const SpmRowElem *spm_elem = rbegin; spm_elem < rend; spm_elem++){
        if (debug)
            std::cerr << "\t" << *spm_elem << "\n";    
        // check if this element contains a pattern
        if (spm_elem->pattern != NULL){
            jmp = this->PreparePat(xs, *spm_elem);
            assert(xs.size() == 0);
            this->AddPattern(*spm_elem, jmp);
            for (long i=0; i < spm_elem->pattern->getSize(); i++){
                this->values[this->values_idx++] = spm_elem->vals[i];
            }
            continue;
        }
        // check if we exceeded the maximum size for a unit
        assert(xs.size() <= CTL_SIZE_MAX);
        if (xs.size() == CTL_SIZE_MAX){
             //std::cerr << "AddXs: max size " << xs.size() << "\n";
             this->AddXs(xs);
        }
        xs.push_back(spm_elem->x);
        this->values[this->values_idx++] = spm_elem->val;
    }
    if (xs.size() > 0){
        //std::cerr << "AddXs: last\n";
        this->AddXs(xs);
    }
}
