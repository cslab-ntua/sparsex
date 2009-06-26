#ifndef CSX_CTL_H__
#define CSX_CTL_H__

extern "C" {
	#include "dynarray.h"
}

namespace csx {

// Ctl Array is a byte-based array storing (compressed)
// index information for the sparse matrix
//
// Format:
//   [flags][size][unit-specific ...]
//
//  flags bits:
//          7  6  5  4  3  2  1  0
//         [ ][ ][ ][ ][ ][ ][ ][ ]
//          |  |  |______________|
//    new row  |        |------------> arbitrary patterns
//          row jmp
//        (if enabled)
#define CTL_NR_BIT 7
#define CTL_RJMP_BIT 6

static inline void set_bit(uint8_t *byte, int bit)
{
	assert(bit >= 0);
	assert(bit < 8);
	*byte |= (1<<bit);
}

// size is 8 bits: 0..255
#define CTL_SIZE_MAX ((1<<8) - 1)

#define CTL_PATTERN_MASK ~(1<<CTL_NR_BIT | 1<<CTL_RJMP_BIT)

class CtlManager
{
public:
	class PatInfo {
	public:
		uint8_t flag; // flags allocated for this pattern
		uint64_t nr;  // number of non-zero elemenets
		PatInfo(uint8_t flag_, uint64_t nr_): flag(flag_), nr(nr_) {}
		PatInfo(): flag(0), nr(0) {}
	};
	typedef std::map<long,PatInfo> PatMap;

	PatMap patterns;
	uint8_t flag_avail; // current available flag
	dynarray_t *ctl_da;
	bool new_row;
	uint64_t last_col;
	SpmIdx *spm;


	CtlManager(SpmIdx *spm_) :flag_avail(0), ctl_da(NULL), spm(spm_) {}

	uint8_t getFlag(long pattern_id, uint64_t nnz);
	uint8_t *mkCtl();
private:
	void doRow(const SpmRowElems &row);
	void AddXs(std::vector<uint64_t> xs);
};

} // end csx namespace

#endif /* CSX_CTL_H__ */
