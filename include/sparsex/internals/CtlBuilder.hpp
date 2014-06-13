/*
 * CtlBuilder.hpp -- Utility class for building the CSX's ctl structure
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CTL_BUILDER_HPP
#define CTL_BUILDER_HPP

#include "sparsex/internals/Allocators.hpp"
#include "sparsex/internals/DynamicArray.hpp"
#include <cstdint>
#include <iostream>

class CtlBuilder
{
public:
    CtlBuilder(size_t size_hint)
        : ctl_da_(CtlDynamicArray(size_hint))
    { }

    size_t GetCtlSize() const
    {
        return ctl_da_.GetSize();
    }

    uint8_t *Finalize()
    {
        // ctl_da_.ShrinkToFit();
        return ctl_da_.TakeElems();
    }

    void AppendCtlHead(bool nr, size_t rowjmp, uint8_t id, uint8_t size,
                       size_t ucol, size_t ucol_size, bool full_colind);
    void AppendVariableInt(unsigned long val);
    void AppendFixedInt(unsigned long val, size_t nr_bytes);

private:
    void AlignCtl(size_t boundary);

#if SPX_USE_NUMA
    typedef DynamicArray<uint8_t, reallocator<uint8_t, NumaAllocator> >
    CtlDynamicArray;
#else
    typedef DynamicArray<uint8_t> CtlDynamicArray;
#endif

    CtlDynamicArray ctl_da_;
};

#endif  // CTL_BUILDER_HPP
