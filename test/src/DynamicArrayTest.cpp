/*
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/*
 * \file DynamicArray.cpp
 * \brief Test dynamic array implementation
 *
 * \author Vasileios Karakasis
 * \date 2013
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/DynamicArray.hpp>
#include <sparsex/internals/Allocators.hpp>
#include <sparsex/internals/Timer.hpp>
#include <iostream>
#include <vector>

using namespace std;
using namespace sparsex::timing;
using namespace sparsex::utilities;

class MyClass
{
public:
    MyClass(size_t a, size_t b)
        : a_(a), b_(b)
    {}

    MyClass()
        : a_(0), b_(0)
    {}

    ~MyClass() {}
private:
    size_t a_, b_;
};

int main(void)
{
    const size_t array_size = 50000000;
//    const size_t array_size = 100;
    const size_t capacity = 10;
    DynamicArray<int, reallocator<int, NumaAllocator> > dynarray(capacity);
    vector<int> vec(capacity);
    Timer timer;

    timer.Start();
    for (size_t i = 0; i < array_size; ++i) {
        //MyClass my(i, i+1);
        dynarray.Append(i);
    }

    timer.Pause();
    cout << "C++ dynamic array impl.: " << timer.ElapsedTime() << " s\n";
    timer.Stop();

    timer.Start();
    for (size_t i = 0; i < array_size; ++i) {
        //MyClass my(i, i+1);
        vec.push_back(i);
    }

    timer.Pause();
    cout << "Vector impl.: " << timer.ElapsedTime() << " s\n";
    timer.Stop();

    StdAllocator &alloc = StdAllocator::GetInstance();
    int *array = new (alloc) int[array_size];

    timer.Start();
    for (size_t i = 0; i < array_size; ++i)
        array[i] = i;

    timer.Pause();
    cout << "Normal array: " << timer.ElapsedTime() << " s\n";
    alloc.Destroy(array, array_size);
    return 0;
}
