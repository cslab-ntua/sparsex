/*
 * Allocators.hpp -- Custom memory allocation
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef ALLOCATORS_HPP
#define ALLOCATORS_HPP

#include "numa_util.h"
#include <cassert>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <new>
#include <vector>

using namespace std;

// not thread safe
class MemoryAllocator
{
public:
    void *Allocate(size_t size, const nothrow_t &)
    {
        return DoAllocate(size);
    }

    void *Allocate(size_t size)
    {
        void *ret = Allocate(size, nothrow);
        if (!ret)
            throw bad_alloc();
        return ret;
    }

    void Deallocate(void *p, size_t size)
    {
        DoDeallocate(p, size);
    }

    void Deallocate(void *p)
    {
        Deallocate(p, 0);
    }

    void *Reallocate(void *p, size_t old_size, size_t new_size,
                     const nothrow_t &)
    {
        return DoReallocate(p, old_size, new_size);
    }

    /**
     *  realloc() contract; see realloc(3)
     */
    void *Reallocate(void *p, size_t old_size, size_t new_size)
    {
        void *ret = Reallocate(p, old_size, new_size, nothrow);
        if (!ret && new_size)
            throw bad_alloc();
        return ret;
    }

    template <typename T>
    void Destroy(T *p)
    {
        if (p) {
            p->~T();
            Deallocate(p, sizeof(T));
        }
    }

    template <typename T>
    void Destroy(T *p, size_t nr_elems)
    {
        if (p) {
            for (size_t i = 0; i < nr_elems; ++i) {
                p[i].~T();
            }

            Deallocate(p, nr_elems*sizeof(T));
        }
    }

protected:
    MemoryAllocator() {};
    virtual ~MemoryAllocator() {};

    // Disable copy/assignment
    MemoryAllocator(const MemoryAllocator &);
    void operator=(const MemoryAllocator &);

private:
    virtual void *DoAllocate(size_t size) throw() = 0;
    virtual void *DoReallocate(void *p, size_t old_size,
                               size_t new_size) throw() = 0;
    virtual void  DoDeallocate(void *p, size_t size) throw() = 0;
};


class StdAllocator
    : public MemoryAllocator
{
public:
    static StdAllocator &GetInstance()
    {
        static StdAllocator instance;
        return instance;
    }

private:
    void *DoAllocate(size_t size) throw()
    {
        return malloc(size);
    }

    void *DoReallocate(void *p, size_t old_size, size_t new_size) throw()
    {
        return realloc(p, new_size);
    }

    void DoDeallocate(void *p, size_t size) throw()
    {
        if (p)
            free(p);
    }

    StdAllocator() {};
    ~StdAllocator() {};

    // Disable copy/assignment
    StdAllocator(const StdAllocator &);
    void operator=(const StdAllocator &);
};

class NumaAllocator
    : public MemoryAllocator
{
public:
    static NumaAllocator &GetInstance()
    {
        static NumaAllocator instance;
        return instance;
    }

    // NUMA-specific allocation methods
    void *AllocateOnNode(size_t size, int node, const nothrow_t &)
    {
        return alloc_onnode(size, node);
    }

    void *AllocateOnNode(size_t size, int node)
    {
        void *ret = AllocateOnNode(size, node, nothrow);
        if (!ret)
            throw bad_alloc();
        return ret;
    }

    void *AllocateInterleaved(size_t size,
                              vector<size_t> parts, vector<int> nodes,
                              const nothrow_t &)
    {
        assert(parts.size() == nodes.size());
        return alloc_interleaved(size, parts.data(), parts.size(),
                                 nodes.data());
    }

    void *AllocateInterleaved(size_t size,
                              vector<size_t> parts, vector<int> nodes)
    {
        void *ret = AllocateInterleaved(size, parts, nodes, nothrow);
        if (!ret)
            throw bad_alloc();
        return ret;
    }
    

private:
    void *DoAllocate(size_t size) throw()
    {
        return numa_alloc_local(size);
    }

    void *DoReallocate(void *p, size_t old_size, size_t new_size) throw()
    {
        if (!new_size)
            DoDeallocate(p, old_size);
        else
            p = numa_realloc(p, old_size, new_size);

        return p;
    }

    void DoDeallocate(void *p, size_t size) throw()
    {
        if (p && size)
            free_interleaved(p, size);
    }

    NumaAllocator() {};
    ~NumaAllocator() {};

    // Disable copy/assignment
    NumaAllocator(const NumaAllocator &);
    void operator=(const NumaAllocator &);
};


// Override new/delete operators using our allocators
void *operator new(size_t size, MemoryAllocator &alloc)
{
    return alloc.Allocate(size);
}

void *operator new[](size_t size, MemoryAllocator &alloc)
{
    return alloc.Allocate(size);
}

void *operator new(size_t size, MemoryAllocator &alloc,
                   const nothrow_t &)
{
    return alloc.Allocate(size, nothrow);
}

void *operator new[](size_t size, MemoryAllocator &alloc,
                     const nothrow_t &)
{
    return alloc.Allocate(size, nothrow);
}

// Provide new operators for NUMA-specific allocations
void *operator new(size_t size, NumaAllocator &alloc, int node)
{
    return alloc.AllocateOnNode(size, node);
}

void *operator new(size_t size, NumaAllocator &alloc, int node,
                   const nothrow_t &)
{
    return alloc.AllocateOnNode(size, node, nothrow);
}

void *operator new[](size_t size, NumaAllocator &alloc, int node)
{
    return alloc.AllocateOnNode(size, node);
}

void *operator new[](size_t size, NumaAllocator &alloc, int node,
                     const nothrow_t &)
{
    return alloc.AllocateOnNode(size, node, nothrow);
}

void *operator new(size_t size, NumaAllocator &alloc,
                   vector<size_t> parts, vector<int> nodes)
{
    return alloc.AllocateInterleaved(size, parts, nodes);
}

void *operator new(size_t size, NumaAllocator &alloc,
                   vector<size_t> parts, vector<int> nodes, const nothrow_t &)
{
    return alloc.AllocateInterleaved(size, parts, nodes, nothrow);
}

void *operator new[](size_t size, NumaAllocator &alloc,
                     vector<size_t> parts, vector<int> nodes)
{
    return alloc.AllocateInterleaved(size, parts, nodes);
}

void *operator new[](size_t size, NumaAllocator &alloc,
                     vector<size_t> parts, vector<int> nodes, const nothrow_t &)
{
    return alloc.AllocateInterleaved(size, parts, nodes, nothrow);
}

void operator delete(void *p, MemoryAllocator &alloc)
{
    return alloc.Deallocate(p);
}

void operator delete[](void *p, MemoryAllocator &alloc)
{
    return alloc.Deallocate(p);
}

template<typename T, typename Alloc = StdAllocator>
class reallocator
{
public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    template<typename U>
    struct rebind {
        typedef reallocator<U> other;
    };

    // All constructors and the destructor are empty, since this class
    // is stateless
    reallocator() throw() {}

    reallocator(const reallocator &other) throw() {}

    template<typename U>
    reallocator(const reallocator<U> &other) throw() {}

    ~reallocator() throw() {}

    pointer address(reference value) const
    {
        return &value;
    }

    const_pointer address(const_reference value) const
    {
        return &value;
    }

    pointer allocate(size_type n, const void *hint = 0)
    {
        return static_cast<pointer>(Alloc::GetInstance().Allocate(n*sizeof(T)));
    }
    
    pointer reallocate(size_type old_n, size_t new_n, void *addr = 0)
    {
        return static_cast<pointer>(
            Alloc::GetInstance().Reallocate(addr, old_n*sizeof(T),
                                            new_n*sizeof(T)));
    }

    pointer reallocate(size_type old_n, size_t new_n, void *addr, 
                       const nothrow_t &) throw()
    {
        return static_cast<pointer>(
            Alloc::GetInstance().Reallocate(addr, old_n*sizeof(T),
                                            new_n*sizeof(T), nothrow));
    }

    void deallocate(pointer p, size_type n)
    {
        Alloc::GetInstance().Deallocate(p, n*sizeof(T));
    }

    size_type max_size() const throw()
    {
        return numeric_limits<size_t>::max() / sizeof(T);
    }

    void construct(pointer p, const_reference val)
    {
        new ((void *) p) value_type(val);
    }

    void destroy(pointer p)
    {
        p->~value_type();
    }

    void destroy(pointer p, size_type n)
    {
        for (size_t i = 0; i < n; ++i)
            p[i].~value_type();
    }
};

template<typename T, typename U>
bool operator==(const reallocator<T> &,
                const reallocator<U> &) throw()
{
    return true;
}

template<typename T, typename U>
bool operator!=(const reallocator<T> &,
                const reallocator<U> &) throw()
{
    return false;
}

#endif  // ALLOCATORS_HPP
