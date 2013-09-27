/* -*- C++ -*-
 *
 * node.h -- Node class.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_NODE_H__
#define CSX_NODE_H__

#include "encodings.h"
#include "sparse_util.h"

#include <iostream>
#include <set>
#include <map>
#include <inttypes.h>

namespace csx {

/**
 *  Keeps data of an encoding sequence.
 */
class Node {
public:
    Node(uint32_t depth);
    ~Node() {}
    
    void PrintNode();

    /**
     *  Ignore the type for the encoding sequence examined.
     *
     *  @param type type which is ignored.
     */
    void Ignore(Encoding::Type type);

    /**
     *  Copies a node to a new one and inserts an extra type in the end of the
     *  encoding sequence.
     *  
     *  @param type   type which is inserted in the end.
     *  @param deltas deltas corresponding to type inserted.
     */
    Node MakeChild(Encoding::Type type, std::set<uint64_t> deltas);

private:
    uint32_t depth_;
    std::map<Encoding::Type, std::set<uint64_t> > deltas_path_;
    Encoding::Type *type_path_;
    Encoding::Type *type_ignore_;
    template<typename IndexType, typename ValueType> friend class EncodingManager;
};

}

#endif  // CSX_NODE_H__
