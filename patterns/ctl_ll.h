#ifndef CTX_CTL_LL_H__
#define CTX_CTL_LL_H__

#include <cassert>

//#include <cstdint>
extern "C" {
	#include <inttypes.h>
}

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

static inline int test_bit(uint8_t *byte, int bit)
{
	assert(bit >= 0);
	assert(bit < 8);
	return (*byte & (1<<bit));
}

// size is 8 bits: 0..255
#define CTL_SIZE_MAX ((1<<8) - 1)

#define CTL_PATTERN_MASK ~(1<<CTL_NR_BIT | 1<<CTL_RJMP_BIT)

#endif // CTL_CTX_LL_H__
