/*
 * node.cc -- Node class implementation.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "node.h"

using namespace csx;

Node::Node(uint32_t depth) : depth_(depth)
{
    uint32_t i;
    
    type_path_ = new IterOrder[depth_];
    type_ignore_ = new IterOrder[XFORM_MAX];
    for (i = 0; i < ((uint32_t) XFORM_MAX); i++)
        type_ignore_[i] = NONE;
}

void Node::PrintNode()
{
    for (uint32_t i = 0; i < depth_; ++i) {
        if (i != 0)
            std::cout << ",";
            
        std::cout << type_path_[i];
    }

    std::cout << std::endl;
    for (uint32_t i = 0; i < depth_; ++i) {
        IterOrder temp_type = type_path_[i];
        std::set<uint64_t>::iterator it = deltas_path_[temp_type].begin();

        if (i != 0)
            std::cout << ",";
            
        std::cout << "{";
        for (uint32_t i = 1;
             i < static_cast<uint32_t>(deltas_path_[temp_type].size());
             ++i) {
            std::cout << *it << ",";
            ++it;
        }

        std::cout << *it << "}";
    }

    std::cout << std::endl;
}

void Node::Ignore(IterOrder type)
{
    uint32_t i = 0;
    
    while (type_ignore_[i] != NONE) {
        assert(type_ignore_[i] != type && "type already ignored");
        ++i;
    }
    
    type_ignore_[i] = type;
}

Node Node::MakeChild(IterOrder type, std::set<uint64_t> deltas)
{
    Node new_node = Node(depth_ + 1);

    for (uint32_t i = 0; i < static_cast<uint32_t>(XFORM_MAX); ++i)
        new_node.type_ignore_[i] = type_ignore_[i];

    for (uint32_t i = 0; i < depth_; ++i) {
        IterOrder temp_type = type_path_[i];
        new_node.type_path_[i] = temp_type;
        new_node.deltas_path_[temp_type] = deltas_path_[temp_type];
    }

    new_node.type_path_[depth_] = type;
    new_node.deltas_path_[type] = deltas;
    return new_node;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
