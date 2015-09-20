/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Element.hpp
 * \brief A generic matrix element of the internal CSX representation
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2015
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_ELEMENT_HPP
#define SPARSEX_INTERNALS_ELEMENT_HPP

#include <sparsex/internals/Encodings.hpp>
#include <sparsex/internals/Utility.hpp>
#include <boost/type_traits.hpp>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <ostream>
#include <string>

using namespace std;
using namespace sparsex::utilities;

namespace sparsex {
namespace csx {

/**
 *  @brief Abstract base class for element marking.
 *
 *  Marking elements can be useful when collecting statistics for the different
 *  encodings of the sparse matrix. For example, you can mark elements in order
 *  to treat them differently when generating or processing statistics. In order
 *  to keep the generic element's implementation as simple as possible, we
 *  introduce the notion of markers. Every element is associated with a marker,
 *  which we can use it to mark it. A marker can accommodate more than a single
 *  mark.
 *
 *  This class provides the basic interface for element markers.
 */
class Marker
{
public:
    /**
     *  Set a mark on this marker.
     *
     *  @param m the mark to set.
     */
    virtual void Mark(unsigned long m) = 0;

    /**
     *  Unset a mark on this marker.
     *
     *  @param m the mark to unset.
     */
    virtual void Unmark(unsigned long m) = 0;

    /**
     *  Clear this marker.
     */ 
    virtual void Unmark() = 0;

    /**
     *  Get the mark of this marker.
     */ 
    virtual unsigned long GetMark() const = 0;

    /**
     *  Check if a specific mark is set on this marker.
     *
     *  @param m the mark to check for.
     */ 
    virtual bool IsMarked(unsigned long m) const = 0;

    /**
     *  Clone this marker.
     */ 
    virtual Marker *Clone() const = 0;

    /**
     *  Destroy this marker.
     */
    virtual ~Marker() {};
};

/**
 *  @brief Mark pattern information
 *
 *  This marker is used for marking the starting element of a pattern and
 *  whether an element will belong to a pattern if the matrix is going to
 *  encoded (see EncodingManager).
 *
 *  Two marks are supported by this marker:
 *  - InPattern: this mark is set for elements that belong to a pattern.
 *  - PatternStart: this mark is set for element at the beginning of a pattern.
 *
 *  @see EncodingManager
 */
class PatternMarker : public Marker
{
public:
    PatternMarker()
        : mark_(0)
    { }

    PatternMarker(const PatternMarker &other)
        : mark_(other.mark_)
    { }

    /**
     *  Set a mark for this marker.
     *  You can use bitwise-OR to set both supported marks, e.g.,
     *  Mark(PatternMarker::InPattern | PatternMarker::PatternStart).
     *
     *  @param m the mark(s) to set on this marker.
     */ 
    void Mark(unsigned long m)
    {
        mark_ |= m;
    }

    /**
     *  Unset a mark for this marker.
     *  You can use bitwise-OR to unset both supported marks, e.g.,
     *  Unmark(PatternMarker::InPattern | PatternMarker::PatternStart).
     *
     *  @param m the mark(s) to unset on this marker.
     */ 
    void Unmark(unsigned long m)
    {
        mark_ &= ~m;
    }

    void Unmark()
    {
        mark_ = 0;
    }

    unsigned long GetMark() const
    {
        return mark_;
    }

    bool IsMarked(unsigned long m) const
    {
        return (mark_ & m) == m;
    }

    PatternMarker *Clone() const
    {
        return new PatternMarker(*this);
    }

    // marker masks
    const static unsigned long InPattern = 1;
    ///< In-pattern marker mask.

    const static unsigned long PatternStart = 2;
    ///< Pattern-start marker mask.

private:
    unsigned long mark_;
};

/**
 *  @brief CSX's generic element.
 *  
 *  This class represents a generic CSX element used in its internal
 *  representation. A generic element is either a single simple non-zero
 *  element, consisting of a row, a column and a value, or a specific CSX
 *  substructure instantiation. In this case, a generic element consists of the
 *  set of non-zero values belonging to the substructure and the coordinates
 *  (row, column) of the upper leftmost element of the instantiation.
 *
 *  This class supplies a move constructor and assignment operator for
 *  efficient copies and assignments from short-lived rvalues.
 *
 *  @see Marker, PatternMarker
 */
template<typename IndexType, typename ValueType>
class Element
{
public:
    /**
     *  Construct a substructure element.
     *
     *  @param row the row of the upper leftmost element in the substructure.
     *  @param col the column of the upper leftmost element in the substructure.
     *  @param vals a forward iterator at the beginning of the substructure's
     *              non-zero values. The values are copied internally.
     *  @param size the size of this instantiation, i.e., the total number of
     *              its comprising elements.
     *  @param inst the encoding instantiation of this generic element.
     *
     *  @throws std::invalid_argument
     *  if size == 0 or vals == NULL.
     *  
     *  @throws std::bad_alloc
     *  if internal values or a default marker cannot be allocated.
     *  
     *  @remark
     *  If size == 1, then this constructor is equivalent to
     *  Element(row, col, *vals).
     *
     *  @remark  
     *  The default marker for this element is PatternMarker.
     */
    template<typename ForwardIterator>
    Element(IndexType row, IndexType col,
            const ForwardIterator &vals, size_t size,
            const Encoding::Instantiation &inst)
        : row_(row),
          col_(col),
          size_(size),
          inst_(inst),
          marker_(new PatternMarker())
    {
        if (!size)
            throw invalid_argument("invalid size");
        check_iterator(vals);
        if (size_ == 1) {
            InitSingleValue(*vals);
        } else {
            vals_ = new ValueType[size];
            copy(vals, vals + size, vals_);
        }
    }

    /**
     *  Construct a single element.
     *
     *  @param row the row of the element.
     *  @param col the column of the element.
     *  @param val the non-zero value of the element.
     *
     *  @throws std::bad_alloc
     *  if default marker cannot be allocated.
     *
     *  @remark
     *  the instantiation of this element is (Encoding::None, 0).
     *  
     *  @remark
     *  The default marker for this element is PatternMarker.
     */
    Element(IndexType row, IndexType col, ValueType val)
        : row_(row),
          col_(col),
          val_(val),
          size_(1),
          inst_(make_pair(Encoding::None, 0)),
          marker_(new PatternMarker())
    { }

    /**
     *  Copy construct an element.
     *
     *  @param other the element to copy from.
     *  
     *  @throws std::bad_alloc
     *  if the internal values or the element's marker cannot be copied.
     */
    Element(const Element<IndexType, ValueType> &other)
        : row_(other.row_),
          col_(other.col_),
          size_(other.size_),
          inst_(other.inst_),
          marker_(other.marker_->Clone())
    {
        if (size_ == 1) {
            val_ = other.val_;
        } else {
            vals_ = new ValueType[size_];
            copy(other.vals_, other.vals_ + size_, vals_);
        }
    }

    /**
     *  Move construct an element.
     *
     *  @param other the element to move from.
     *  @throws This function offers no-throw guarantee.
     *  @remark
     *  The internal values and marker are stolen from other. At exit other's
     *  values and marker are set to 0.
     */
    Element(Element<IndexType, ValueType> &&other)
        : row_(other.row_),
          col_(other.col_),
          size_(other.size_),
          inst_(other.inst_),
          marker_(other.marker_)    // steal the marker from the rvalue
    {
        if (size_ == 1) {
            val_ = other.val_;
        } else {
            vals_ = other.vals_;
            other.vals_ = 0;
        }

        other.marker_ = 0;
    }

    /**
     *  Destroy this element.
     */
    ~Element()
    {
        DestroyVals();
        DestroyMarker();
    }

    /**
     *  Get the size of this element.
     *
     *  @return the size of this element.
     */ 
    size_t GetSize() const
    {
        return size_;
    }

    /**
     *  Get the row of this element.
     *
     *  @return the row of this element.
     */ 
    IndexType GetRow() const
    {
        return row_;
    }

    /**
     *  Get the column of this element.
     *
     *  @return the column of this element.
     */ 
    IndexType GetCol() const
    {
        return col_;
    }

    /**
     *  Get the coordinates of this element as a pair.
     *
     *  This call is equivalent to std::make_pair(GetRow(), GetCol())
     *
     *  @return the coordinates of this element.
     */ 
    pair<IndexType, IndexType> GetCoordinates() const
    {
        return make_pair(row_, col_);
    }

    /**
     *  Checks whether this element is a substructure.
     *
     *  @return true if this element is a substructure, false otherwise.
     */
    bool IsPattern() const
    {
        // Only the instantiation can show whether we are a stray element or not
        // if inst_.first == Encoding::None, we may be in a delta unit
        return inst_.second != 0;
    }

    /**
     *  Get the type of this element
     *
     *  @return the type of this element. If this is a single element, then
     *  Encoding::None is returned.
     */
    Encoding::Type GetType() const
    {
        return inst_.first;
    }

    const Encoding::Instantiation &GetInstantiation() const
    {
        return inst_;
    }

    /**
     *  Get the values of this substructure element.
     *
     *  @return a const reference to this element's values.
     *
     *  @throws std::invalid_argument if this is not a substructure element.
     *
     *  @remark  
     *  This method can only be invoked on substructure elements.
     */
    const ValueType &GetValues() const
    {
        if (size_ == 1)
            throw invalid_argument("cannot GetValues() on single element");
        return *vals_;
    }

    /**
     *  Get the value of this single element.
     *
     *  @return the value of this element.
     *  @throws This function offers no-throw guarantee.
     *  @remark
     *  The returned reads correctly on single elements.
     */
    ValueType GetValue() const
    {
        return val_;
    }

    /**
     *  Get the marker of this element
     *
     *  @return a reference to the marker of this element.
     */
    Marker &GetMarker() const
    {
        return *marker_;
    }

    /**
     *  Transform the coordinates of this element.
     */
    void Transform(pair<IndexType, IndexType> new_coord)
    {
        row_ = new_coord.first;
        col_ = new_coord.second;
    }

    /**
     *  Set the marker of this element.
     *
     *  The supplied marker m is copied internally.
     *
     *  @throw std::bad_alloc if copying of m failed.
     */
    void SetMarker(const Marker &m)
    {
        DestroyMarker();
        marker_ = m.Clone();
    }

    // Relational operators: equality is checked only lexicographically

    /**
     *  Equality operator.
     *
     *  Equality is checked only on coordinates, i.e., *this == rhs, if
     *  this->GetRow() == rhs.GetRow() && this->GetCol() == rhs.GetCol().
     *
     *  @return true if this element is equal to rhs, false otherwise.
     */
    bool operator==(const Element<IndexType, ValueType> &rhs) const
    {
        return row_ == rhs.row_ && col_ == rhs.col_;
    }

    bool operator!=(const Element<IndexType, ValueType> &rhs) const
    {
        return !(*this == rhs);
    }

    /**
     *  Less than operator.
     *
     *  The elements are compared lexicographically based on their coordinates,
     *  i.e., *this < rhs, if this->GetRow() < rhs.GetRow() || (this->GetRow()
     *  == rhs.GetRow() && this->GetCol() < rhs.GetCol()).
     *
     *  @return true if this is less than rhs, false otherwise.
     */ 
    bool operator<(const Element<IndexType, ValueType> &rhs) const
    {
        return row_ < rhs.row_ || (row_ == rhs.row_ && col_ < rhs.col_);
    }

    bool operator<=(const Element<IndexType, ValueType> &rhs) const
    {
        return *this < rhs || *this == rhs;
    }

    bool operator>(const Element<IndexType, ValueType> &rhs) const
    {
        return !(*this <= rhs);
    }

    bool operator>=(const Element<IndexType, ValueType> &rhs) const
    {
        return *this > rhs || *this == rhs;
    }

    template<typename I, typename V>
    friend ostream &operator<<(ostream &, const Element<I, V> &);
    template<typename I, typename V>
    friend void swap(Element<I, V> &, Element<I, V> &);

    /**
     *  Assignment operator using copy and swap idiom.
     *
     *  @return this updated instance
     *
     *  @throws std::bad_alloc if other cannot be copy constructed locally.
     */
    // Element<IndexType, ValueType> &operator=(
    //     Element<IndexType, ValueType> other)
    // {
    //     swap(*this, other);
    //     return *this;
    // }

    /**
     *  Copy assignment operator.
     *
     *  @return this updated instance
     *
     *  @throws std::bad_alloc if other cannot be copy constructed locally.
     */
    Element<IndexType, ValueType> &operator=(
        const Element<IndexType, ValueType> &other)
    {
        row_ = other.row_;
        col_ = other.col_;
        size_ = other.size_;
        inst_ = other.inst_;
        marker_ = other.marker_->Clone();
        if (size_ == 1) {
            val_ = other.val_;
        } else {
            vals_ = new ValueType[size_];
            copy(other.vals_, other.vals_ + size_, vals_);
        }

        return *this;
    }

    /**
     *  Move assignment operator.
     *
     *  @return this updated instance
     */
    Element<IndexType, ValueType> &operator=(
        Element<IndexType, ValueType>&& other)
    {
        row_ = other.row_;
        col_ = other.col_;
        size_ = other.size_;
        inst_ = other.inst_;
        marker_ = other.marker_;
        if (size_ == 1) {
            val_ = other.val_;
        } else {
            vals_ = other.vals_;
            other.vals_ = 0;
        }

        other.marker_ = 0;
        return *this;
    }

private:
    void InitSingleValue(ValueType val)
    {
        size_ = 1;
        val_ = val;
        inst_ = make_pair(Encoding::None, 0);
    }

    void DestroyVals()
    {
        if (size_ > 1 && vals_)
            delete[] vals_;
    }
    
    void DestroyMarker()
    {
        if (marker_)
            delete marker_;
    }

    IndexType row_;
    IndexType col_;

    // Use a union for either the value of single element or the values
    // contained in pattern, in order to avoid unnecessary allocations and
    // deep copies for single elements
    union {
        ValueType val_;
        ValueType *vals_;
    };

    size_t size_;
    Encoding::Instantiation inst_;
    Marker *marker_;    // element marker
};

template<typename IndexType, typename ValueType>
ostream &operator<<(ostream &out, const Element<IndexType, ValueType> &e)
{
    if (e.size_ == 1) {
        out << "(" << e.row_ << ", " << e.col_ << ", " << e.val_ << ")";
    } else {
        out << "(" << e.row_ << ", " << e.col_ << ", {";
        for (size_t i = 0; i < e.size_ - 1; ++i)
            out << e.vals_[i] << ", ";
        Encoding enc(e.inst_.first);
        out << e.vals_[e.size_-1] << "}) [s: " << e.size_
            << ", inst: (" << enc << ", " << e.inst_.second << ")]";
    }

    return out;
}

template<typename IndexType, typename ValueType>
void swap(Element<IndexType, ValueType> &a, Element<IndexType, ValueType> &b)
{
    using std::swap;
    swap(a.row_, b.row_);
    swap(a.col_, b.col_);
    if (sizeof(a.vals_) >= sizeof(a.val_))
        swap(a.vals_, b.vals_);
    else
        swap(a.val_, b.val_);

    swap(a.size_, b.size_);
    swap(a.inst_, b.inst_);
    swap(a.marker_, b.marker_);
}

/**
 *  Get the last column of elem's non-zero values based on a specific iteration
 *  order.
 *
 *  If the desired iteration order is different from elem's, then the
 *  elem.GetCol() is returned. Otherwise, the column of the last non-zero value
 *  in the elem's substructure is returned.
 *
 *  @param elem the element whose last column will be computed
 *  @param order the iteration order according to which the last column will
 *               be computed
 *  @return the last column index of this generic element.
 */
template<typename IndexType, typename ValueType>
size_t GetLastCol(const Element<IndexType, ValueType> &elem,
                  Encoding::Type order)
{
    size_t ret = elem.GetCol();
    if (order == elem.GetType())
        ret += (elem.GetSize() - 1) * elem.GetInstantiation().second;

    return ret;
}

} // end of namespace csx
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_ELEMENT_HPP
