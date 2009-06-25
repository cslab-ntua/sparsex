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

// size is 8 bits: 0..255
#define CTL_SIZE_MAX ((1<<8) - 1)

#define CTL_PATTERN_MASK ~(1<<CTL_NR_BIT | 1<<CTL_RJMP_BIT)

class CtlManager
{
public:
	typedef class {
	public:
		uint8_t flag; // flags allocated for this pattern
		uint64_t nr;  // number of non-zero elemenets
		PatInfo(uint8_t flag_, uint64_t nr_): flag(flag_), nr(nr_) {}
	} PatInfo;
	typedef map<long,PatInfo> PatMap;

	PatMap patterns;
	uint8_t flag_avail; // current available flag
	dynarray_t *ctl_da;
	SpmIdx *spm;


	CtlManager(SpmIdx *spm_) :flag_avail(0), ctl_da(NULL), spm(spm_)
	{
		this->ctl_da = dynarray_create(sizeof(uint8_t), 512)
	}

	getFlag(long pattern_id);
	uint8_t *mkCtl();
};

} // end csx namespace

#endif /* CSX_CTL_H__ */
