/*
 * CtlBuilder.cpp -- Utility class for building the CSX's ctl structure
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "sparsex/internals/CtlBuilder.hpp"
#include "sparsex/internals/CtlUtil.hpp"

static inline void SetBit(uint8_t *byte, int bit)
{
    assert(bit >= 0);
    assert(bit < 8);
    *byte |= (1<<bit);
}

void CtlBuilder::AppendVariableInt(unsigned long val)
{
    const unsigned shift = 7;
    uint8_t byte = 0;

    for (;;) {
        byte = (val & ((1 << shift) - 1));
        if (val < (1 << shift)) {
            ctl_da_.Append(byte);
            break;
        }

        byte |= (1 << shift);
        ctl_da_.Append(byte);
        val >>= shift;
    }
}

void CtlBuilder::AppendFixedInt(unsigned long val, size_t nr_bytes)
{
    assert(nr_bytes <= sizeof(val));
#ifdef PTR_ALIGN
    AlignCtl(nr_bytes);
#endif
    for (size_t i = 0; i < nr_bytes; ++i) {
        uint8_t byte = *((uint8_t *) &val + i);
        ctl_da_.Append(byte);
    }
}

void CtlBuilder::AppendCtlHead(bool nr, size_t rowjmp, uint8_t id, uint8_t size,
                               size_t ucol, size_t ucol_size, bool full_colind)
{
    uint8_t flag = id;
    if (nr)
        SetBit(&flag, CTL_NR_BIT);
    if (rowjmp)
        SetBit(&flag, CTL_RJMP_BIT);

    ctl_da_.Append(flag);
    ctl_da_.Append(size);
    if (rowjmp)
        AppendVariableInt(rowjmp);

    if (full_colind)
        AppendFixedInt(ucol, ucol_size);
    else
        AppendVariableInt(ucol);
}

void CtlBuilder::AlignCtl(size_t boundary)
{
    // Explicit padding must be added to ensure right pointer alignment
    // when ctl is dereferenced by the SpMV routine
    const uint8_t *elems = ctl_da_.GetElems();
    const uint8_t *next = &elems[ctl_da_.GetSize()];
	size_t nr_padd = (boundary - ((size_t) next & (boundary-1))) & (boundary-1);
    for (size_t i = 0; i < nr_padd; ++i)
        ctl_da_.Append(0);
}
