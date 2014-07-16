/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Utility.hpp
 * \brief Several generic utility functions
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_UTILITY_HPP
#define SPARSEX_INTERNALS_UTILITY_HPP

#include <boost/type_traits.hpp>

using namespace std;

namespace sparsex {
namespace utilities {

template<bool property>
struct math_impl
{
    template<typename T>
    static T do_gcd(T a, T b);

    template <typename T>
    static T do_abs(T a);

    template <typename T>
    static T do_iceil(T a, T b);
};

template<typename T>
T gcd(T a, T b)
{
    return math_impl<boost::is_integral<T>::value>::do_gcd(a, b);
}

// Integer ceiling
template<typename T>
T iceil(T a, T b)
{
    return math_impl<boost::is_integral<T>::value>::do_iceil(a, b);
}

template<typename T>
T abs(T a)
{
    return math_impl<boost::is_arithmetic<T>::value>::do_abs(a);
}

template<typename T>
T lcm(T a, T b)
{
    return abs(a*b) / gcd(a, b);
}

template<>
struct math_impl<true> 
{
    template<typename T>
    static T do_gcd(T a, T b)
    {
        while (b) {
            T t = b;
            b = a % b;
            a = t;
        }

        return a;
    }

    template<typename T>
    static T do_iceil(T a, T b)
    {
        return a / b + (a % b != 0);
    }

    template<typename T>
    static T do_abs(T a)
    {
        return (a < 0) ? -a : a;
    }
};

template<bool is_pointer>
struct iterator_checker
{
    template<typename Iterator>
    static void do_check(const Iterator &) {}
};

template<typename Iterator>
void check_iterator(const Iterator &i)
{
    iterator_checker<boost::is_pointer<Iterator>::value>::do_check(i);
}

template<>
struct iterator_checker<true>
{
    template<typename Iterator>
    static void do_check(const Iterator &i)
    {
        if (!i)
            throw invalid_argument("invalid iterator");
    }
};

} // end of namespace utilities
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_UTILITY_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
