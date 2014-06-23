/*
 * Node.cpp -- Node class implementation.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "sparsex/internals/Node.hpp"

using namespace csx;
using namespace std;

Node::Node(uint32_t depth)
  : depth_(depth)
{
    type_path_ = new Encoding::Type[depth_];
    type_ignore_ = new Encoding::Type[Encoding::Max];
    for (Encoding::Type t = Encoding::None; t < Encoding::Max; t++)
        type_ignore_[t] = Encoding::None;
}

void Node::PrintNode()
{
    for (uint32_t i = 0; i < depth_; ++i) {
        if (i != 0)
            cout << ",";
            
        cout << type_path_[i];
    }

    cout << endl;
    for (uint32_t i = 0; i < depth_; ++i) {
        Encoding::Type temp_type = type_path_[i];
        set<uint64_t>::iterator it = deltas_path_[temp_type].begin();

        if (i != 0)
            cout << ",";
            
        cout << "{";
        for (uint32_t i = 1;
             i < static_cast<uint32_t>(deltas_path_[temp_type].size());
             ++i) {
            cout << *it << ",";
            ++it;
        }

        cout << *it << "}";
    }

    cout << endl;
}

void Node::Ignore(Encoding::Type type)
{
    uint32_t i = 0;
    while (type_ignore_[i] != Encoding::None) {
        assert(type_ignore_[i] != type && "type already ignored");
        ++i;
    }
    
    type_ignore_[i] = type;
}

Node Node::MakeChild(Encoding::Type type, set<uint64_t> deltas)
{
    Node new_node = Node(depth_ + 1);

    for (Encoding::Type t = Encoding::None; t < Encoding::Max; ++t)
        new_node.type_ignore_[t] = type_ignore_[t];

    for (size_t i = 0; i < depth_; ++i) {
        Encoding::Type temp_type = type_path_[i];
        new_node.type_path_[i] = temp_type;
        new_node.deltas_path_[temp_type] = deltas_path_[temp_type];
    }

    new_node.type_path_[depth_] = type;
    new_node.deltas_path_[type] = deltas;
    return new_node;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
