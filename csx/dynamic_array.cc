/* -*- C++ -*-
 *
 * dynamicarray.h -- NUMA aware allocator
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include "dynamic_array.h"
#include "allocators.h"
#include "timer.h"

#include <iostream>
#include <vector>

using namespace std;

class MyClass
{
public:
    MyClass(std::size_t a, std::size_t b)
        : a_(a), b_(b)
    {}

    MyClass()
        : a_(0), b_(0)
    {}

    ~MyClass() {}
private:
    std::size_t a_, b_;
};

int main(void)
{
    const size_t array_size = 100000000;
//    const size_t array_size = 100;
    const size_t capacity = 10;
    DynamicArray<MyClass> dynarray(capacity);
    vector<MyClass> vec(capacity);

    xtimer_t timer;

    timer_init(&timer);
    timer_start(&timer);
    for (size_t i = 0; i < array_size; ++i) {
        MyClass my(i, i+1);
        dynarray.Append(my);
    }
    timer_pause(&timer);
    std::cout << "Dynamic array impl.: " << timer_secs(&timer) << " s\n";

    timer_init(&timer);
    timer_start(&timer);
    for (size_t i = 0; i < array_size; ++i) {
        MyClass my(i, i+1);
        vec.push_back(my);
    }

    timer_pause(&timer);
    std::cout << "Vector impl.: " << timer_secs(&timer) << " s\n";

    StdAllocator &alloc = StdAllocator::GetInstance();
    int *array = new (alloc) int[array_size];

    timer_init(&timer);
    timer_start(&timer);
    for (size_t i = 0; i < array_size; ++i)
        array[i] = i;
    timer_pause(&timer);
    std::cout << "Normal array: " << timer_secs(&timer) << " s\n";
    alloc.Destroy(array, array_size);
    return 0;
}
