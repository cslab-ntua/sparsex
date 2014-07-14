/*
 * \file CsrIterator.hpp
 *
 * \brief An iterator that simultaneously iterates through two arrays
 * 
 * Based on iterator_facade by boost, compliant with std::sort.
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_CSR_ITERATOR_HPP
#define SPARSEX_INTERNALS_CSR_ITERATOR_HPP

#include <iostream>
#include <algorithm>
#include <inttypes.h>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/variant/get.hpp>
#include <boost/tuple/tuple.hpp>

using namespace std;

namespace sparsex {
namespace io {

template <typename SortIter, typename PermIter> 
struct CSR_iterator_traits
{ 
    typedef boost::tuple< 
        typename std::iterator_traits<SortIter>::value_type, 
        typename std::iterator_traits<PermIter>::value_type > 
        value_type; 

    typedef boost::tuple< 
        typename std::iterator_traits<SortIter>::value_type&, 
        typename std::iterator_traits<PermIter>::value_type& > 
        ref_type; 
};

template<typename SortIter, typename PermIter>
class CSR_iterator
    : public boost::iterator_facade<
    CSR_iterator<SortIter, PermIter>,
    typename CSR_iterator_traits<SortIter, PermIter>::value_type,
    std::random_access_iterator_tag,
    typename CSR_iterator_traits<SortIter, PermIter>::ref_type,
    typename std::iterator_traits<SortIter>::difference_type
    > 
{
public:

    CSR_iterator() {}
    CSR_iterator(SortIter si, PermIter pi) : si_(si), pi_(pi) {}

    SortIter si_;
    PermIter pi_;

private:

    friend class boost::iterator_core_access;

    typename CSR_iterator_traits<SortIter, PermIter>::ref_type dereference() const
    {
        return (typename CSR_iterator_traits<SortIter, PermIter>
                ::ref_type(*si_, *pi_)); 
    }

    bool equal(CSR_iterator const& other) const
    { 
        return (si_ == other.si_); 
    }

    void increment()
    {
        ++si_;
        ++pi_;
    }

    void decrement()
    {
        --si_;
        --pi_;
    }

    void advance(typename std::iterator_traits<SortIter>::difference_type n)
    {
       si_ += n;
       pi_ += n;
    }

    typename std::iterator_traits<SortIter>::difference_type distance_to(CSR_iterator const& other) const
    {
        return (other.si_ - si_); 
    }
};

template <typename SortIter, typename PermIter> 
CSR_iterator<SortIter, PermIter>
get_CSR_iterator(SortIter si, PermIter pi) 
{ 
    return CSR_iterator<SortIter, PermIter>(si, pi); 
}

template <typename SortIter, typename PermIter> 
struct CSR_Comp : public std::binary_function<
    typename CSR_iterator_traits<SortIter, PermIter>::value_type,
    typename CSR_iterator_traits<SortIter, PermIter>::value_type, bool> 
{
    typedef typename CSR_iterator_traits<SortIter, PermIter>::value_type T; 
    
    bool operator()(const T& t1, const T& t2) 
    { 
        return (boost::get<0>(t1) < boost::get<0>(t2)); 
    } 
};

} // end of namespace io
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_CSR_ITERATOR_HPP
