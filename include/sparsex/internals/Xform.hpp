/*
 * \file Xform.hpp
 *
 * \brief Matrix coordinate transformations
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_XFORM_HPP
#define SPARSEX_INTERNALS_XFORM_HPP

#include <sparsex/internals/Element.hpp>
#include <sparsex/internals/Encodings.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/function.hpp>
#include <cassert>

namespace bll = boost::lambda;
using namespace std;

namespace sparsex {
namespace csx {

// Coordinate transformation functions
// The variant with the pair argument is used for chaining transformations

template<typename IndexType>
pair<IndexType, IndexType> XformHoriz(IndexType r, IndexType c,
                                      IndexType nr_rows = 0,
                                      IndexType nr_cols = 0)
{
    // Identinty function
    return make_pair(r, c);
}

template<typename IndexType>
pair<IndexType, IndexType> XformHoriz(pair<IndexType, IndexType> coord,
                                      IndexType nr_rows = 0,
                                      IndexType nr_cols = 0)
{
    return XformHoriz<IndexType>(coord.first, coord.second);
}


template<typename IndexType>
pair<IndexType, IndexType> RevXformHoriz(IndexType r, IndexType c,
                                         IndexType nr_rows = 0,
                                         IndexType nr_cols = 0)
{
    // Identinty transformation
    return make_pair(r, c);
}


template<typename IndexType>
pair<IndexType, IndexType> RevXformHoriz(pair<IndexType, IndexType> coord,
                                         IndexType nr_rows = 0,
                                         IndexType nr_cols = 0)
{
    return RevXformHoriz<IndexType>(coord.first, coord.second);
}

template<typename IndexType>
pair<IndexType, IndexType> XformVert(IndexType r, IndexType c,
                                     IndexType nr_rows = 0,
                                     IndexType nr_cols = 0)
{
    return make_pair(c, r);
}

template<typename IndexType>
pair<IndexType, IndexType> XformVert(pair<IndexType, IndexType> coord,
                                     IndexType nr_rows = 0,
                                     IndexType nr_cols = 0)
{
    return XformVert<IndexType>(coord.first, coord.second);
}

template<typename IndexType>
pair<IndexType, IndexType> RevXformVert(IndexType r, IndexType c,
                                        IndexType nr_rows = 0,
                                        IndexType nr_cols = 0)
{
    return XformVert<IndexType>(r, c);
}

template<typename IndexType>
pair<IndexType, IndexType> RevXformVert(pair<IndexType, IndexType> coord,
                                        IndexType nr_rows = 0,
                                        IndexType nr_cols = 0)
{
    return RevXformVert<IndexType>(coord.first, coord.second);
}

template<typename IndexType>
pair<IndexType, IndexType> XformDiag(IndexType r, IndexType c,
                                     IndexType nr_rows,
                                     IndexType nr_cols = 0)
{
    assert(nr_rows + c - r > 0);
    return make_pair(nr_rows + c - r, (c < r) ? c : r);
}

template<typename IndexType>
pair<IndexType, IndexType> XformDiag(pair<IndexType, IndexType> coord,
                                     IndexType nr_rows,
                                     IndexType nr_cols = 0)
{
    return XformDiag<IndexType>(coord.first, coord.second, nr_rows);
}

template<typename IndexType>
pair<IndexType, IndexType> RevXformDiag(IndexType r, IndexType c,
                                        IndexType nr_rows,
                                        IndexType nr_cols = 0)
{
    if (r < nr_rows)
        return make_pair(nr_rows + c - r, c);
    else
        return make_pair(c, r + c - nr_rows);
}

template<typename IndexType>
pair<IndexType, IndexType> RevXformDiag(pair<IndexType, IndexType> coord,
                                        IndexType nr_rows,
                                        IndexType nr_cols = 0)
{
    return RevXformDiag<IndexType>(coord.first, coord.second, nr_rows);
}

template<typename IndexType>
pair<IndexType, IndexType> XformAntiDiag(IndexType r, IndexType c,
                                         IndexType nr_rows,
                                         IndexType nr_cols)
{
    IndexType new_r = r + c - 1;
    return make_pair(new_r, (new_r <= nr_cols) ? r : nr_cols - c + 1);
}

template<typename IndexType>
pair<IndexType, IndexType> XformAntiDiag(pair<IndexType, IndexType> coord,
                                         IndexType nr_rows, IndexType nr_cols)
{
    return XformAntiDiag<IndexType>(coord.first, coord.second,
                                    nr_rows, nr_cols);
}

template<typename IndexType>
pair<IndexType, IndexType> RevXformAntiDiag(IndexType r, IndexType c,
                                            IndexType nr_rows,
                                            IndexType nr_cols)
{
    if (r <= nr_cols)
        return make_pair(c, r - c + 1);
    else
        return make_pair(r + c - nr_cols, nr_cols - c + 1);
}

template<typename IndexType>
pair<IndexType, IndexType> RevXformAntiDiag(pair<IndexType, IndexType> coord,
                                            IndexType nr_rows,
                                            IndexType nr_cols)
{
    return RevXformAntiDiag<IndexType>(coord.first, coord.second,
                                       nr_rows, nr_cols);
}

template<typename IndexType, IndexType R>
pair<IndexType, IndexType> XformBlockRow(IndexType r, IndexType c,
                                         IndexType nr_rows = 0,
                                         IndexType nr_cols = 0)
{
    return make_pair((r - 1) / R + 1,
                     (r - 1) % R + R*(c - 1) + 1);
}

template<typename IndexType, IndexType R>
pair<IndexType, IndexType> XformBlockRow(pair<IndexType, IndexType> coord,
                                         IndexType nr_rows = 0,
                                         IndexType nr_cols = 0)
{
    return XformBlockRow<IndexType, R>(coord.first, coord.second);
}

template<typename IndexType, IndexType R>
pair<IndexType, IndexType> RevXformBlockRow(IndexType r, IndexType c,
                                            IndexType nr_rows = 0,
                                            IndexType nr_cols = 0)
{
    return make_pair(R*(r - 1) + (c - 1) % R + 1,
                     (c - 1) / R + 1);
}

template<typename IndexType, IndexType R>
pair<IndexType, IndexType> RevXformBlockRow(pair<IndexType, IndexType> coord,
                                            IndexType nr_rows = 0,
                                            IndexType nr_cols = 0)
{   
    return RevXformBlockRow<IndexType, R>(coord.first, coord.second);
}

template<typename IndexType, IndexType C>
pair<IndexType, IndexType> XformBlockCol(IndexType r, IndexType c,
                                         IndexType nr_rows = 0,
                                         IndexType nr_cols = 0)
{
    pair<IndexType, IndexType> vert = XformVert<IndexType>(r, c);
    return XformBlockRow<IndexType, C>(vert.first, vert.second);
}

template<typename IndexType, IndexType C>
pair<IndexType, IndexType> XformBlockCol(pair<IndexType, IndexType> coord,
                                         IndexType nr_rows = 0,
                                         IndexType nr_cols = 0)
{
    return XformBlockCol<IndexType, C>(coord.first, coord.second);
}

template<typename IndexType, IndexType C>
pair<IndexType, IndexType> RevXformBlockCol(IndexType r, IndexType c,
                                            IndexType nr_rows = 0,
                                            IndexType nr_cols = 0)
{
    pair<IndexType, IndexType> rblock = RevXformBlockRow<IndexType, C>(r, c);
    return RevXformVert<IndexType>(rblock.first, rblock.second);
}

template<typename IndexType, IndexType C>
pair<IndexType, IndexType> RevXformBlockCol(pair<IndexType, IndexType> coord,
                                            IndexType nr_rows = 0,
                                            IndexType nr_cols = 0)
{
    return RevXformBlockCol<IndexType, C>(coord.first, coord.second);
}

// Transform function type
template<typename IndexType>
struct TransformFn {
    typedef boost::function<
        pair<IndexType, IndexType>
        (pair<IndexType, IndexType>, IndexType, IndexType)> type;
};

template<typename IndexType>
typename TransformFn<IndexType>::type GetXformFnFromHoriz(Encoding::Type to)
{
    typename TransformFn<IndexType>::type ret;

    switch (to) {
    case Encoding::Horizontal:
        ret = bll::bind(XformHoriz<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::Vertical:
        ret = bll::bind(XformVert<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::Diagonal:
        ret = bll::bind(XformDiag<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::AntiDiagonal:
        ret = bll::bind(XformAntiDiag<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow1:
        ret = bll::bind(XformBlockRow<IndexType, 1>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow2:
        ret = bll::bind(XformBlockRow<IndexType, 2>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow3:
        ret = bll::bind(XformBlockRow<IndexType, 3>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow4:
        ret = bll::bind(XformBlockRow<IndexType, 4>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow5:
        ret = bll::bind(XformBlockRow<IndexType, 5>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow6:
        ret = bll::bind(XformBlockRow<IndexType, 6>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow7:
        ret = bll::bind(XformBlockRow<IndexType, 7>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow8:
        ret = bll::bind(XformBlockRow<IndexType, 8>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol1:
        ret = bll::bind(XformBlockCol<IndexType, 1>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol2:
        ret = bll::bind(XformBlockCol<IndexType, 2>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol3:
        ret = bll::bind(XformBlockCol<IndexType, 3>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol4:
        ret = bll::bind(XformBlockCol<IndexType, 4>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol5:
        ret = bll::bind(XformBlockCol<IndexType, 5>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol6:
        ret = bll::bind(XformBlockCol<IndexType, 6>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol7:
        ret = bll::bind(XformBlockCol<IndexType, 7>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol8:
        ret = bll::bind(XformBlockCol<IndexType, 8>, bll::_1, bll::_2, bll::_3);
        break;
    default:
        assert(false && "[BUG] `to' is not handled");
    }

    return ret;
}

template<typename IndexType>
typename TransformFn<IndexType>::type GetXformFnToHoriz(Encoding::Type from)
{
    typename TransformFn<IndexType>::type ret;

    switch (from) {
    case Encoding::Horizontal:
        ret = bll::bind(RevXformHoriz<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::Vertical:
        ret = bll::bind(RevXformVert<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::Diagonal:
        ret = bll::bind(RevXformDiag<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::AntiDiagonal:
        ret = bll::bind(RevXformAntiDiag<IndexType>, bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow1:
        ret = bll::bind(RevXformBlockRow<IndexType, 1>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow2:
        ret = bll::bind(RevXformBlockRow<IndexType, 2>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow3:
        ret = bll::bind(RevXformBlockRow<IndexType, 3>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow4:
        ret = bll::bind(RevXformBlockRow<IndexType, 4>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow5:
        ret = bll::bind(RevXformBlockRow<IndexType, 5>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow6:
        ret = bll::bind(RevXformBlockRow<IndexType, 6>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow7:
        ret = bll::bind(RevXformBlockRow<IndexType, 7>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockRow8:
        ret = bll::bind(RevXformBlockRow<IndexType, 8>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol1:
        ret = bll::bind(RevXformBlockCol<IndexType, 1>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol2:
        ret = bll::bind(RevXformBlockCol<IndexType, 2>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol3:
        ret = bll::bind(RevXformBlockCol<IndexType, 3>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol4:
        ret = bll::bind(RevXformBlockCol<IndexType, 4>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol5:
        ret = bll::bind(RevXformBlockCol<IndexType, 5>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol6:
        ret = bll::bind(RevXformBlockCol<IndexType, 6>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol7:
        ret = bll::bind(RevXformBlockCol<IndexType, 7>,
                        bll::_1, bll::_2, bll::_3);
        break;
    case Encoding::BlockCol8:
        ret = bll::bind(RevXformBlockCol<IndexType, 8>,
                        bll::_1, bll::_2, bll::_3);
        break;
    default:
        assert(false && "[BUG] `from' is not handled");
    }

    return ret;
}

template<typename IndexType>
typename TransformFn<IndexType>::type GetXformFn(Encoding::Type from,
                                                 Encoding::Type to)
{
    typename TransformFn<IndexType>::type ret;

    if (from == to)
        ret = bll::bind(XformHoriz<IndexType>, bll::_1, bll::_2, bll::_3);
    else if (from == Encoding::Horizontal)
        ret = GetXformFnFromHoriz<IndexType>(to);
    else if (to == Encoding::Horizontal)
        ret = GetXformFnToHoriz<IndexType>(from);
    else {
        typename TransformFn<IndexType>::type from_horiz =
            GetXformFnFromHoriz<IndexType>(to);
        typename TransformFn<IndexType>::type to_horiz =
            GetXformFnToHoriz<IndexType>(from);
        ret = bll::bind(from_horiz,
                        bll::bind(to_horiz, bll::_1, bll::_2, bll::_3),
                        bll::_2, bll::_3);
    }

    return ret;
}

}   // end of namespace csx
}   // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_XFORM_HPP
