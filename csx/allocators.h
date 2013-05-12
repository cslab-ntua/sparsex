/* -*- C++ -*-
 *
 * allocators.h -- Custom memory allocation
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef ALLOC_H
#define ALLOC_H

#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <new>

// not thread safe
class MemoryAllocator
{
public:
    void *Allocate(std::size_t size, const std::nothrow_t &)
    {
        return DoAllocate(size);
    }

    void *Allocate(std::size_t size)
    {
        void *ret = Allocate(size, std::nothrow);
        if (!ret)
            throw std::bad_alloc();
        return ret;
    }

    void Deallocate(void *p, std::size_t size)
    {
        DoDeallocate(p, size);
    }

    void Deallocate(void *p)
    {
        Deallocate(p, 0);
    }

    void *Reallocate(void *p, std::size_t size, const std::nothrow_t &)
    {
        return DoReallocate(p, size);
    }

    void *Reallocate(void *p, std::size_t size)
    {
        void *ret = Reallocate(p, size, std::nothrow);
        if (!ret)
            throw std::bad_alloc();
        return ret;
    }

    template <typename T>
    void Destroy(T *p)
    {
        if (p) {
            p->~T();
            std::cout << "Size: " << sizeof(T) << "\n";
            Deallocate(p, sizeof(T));
        }
    }

    template <typename T>
    void Destroy(T *p, std::size_t nr_elems)
    {
        if (p) {
            for (std::size_t i = 0; i < nr_elems; ++i) {
                p->~T();
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
    virtual void *DoAllocate(std::size_t size) throw() = 0;
    virtual void *DoReallocate(void *p, std::size_t size) throw() = 0;
    virtual void  DoDeallocate(void *p, std::size_t size) throw() = 0;
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
    void *DoAllocate(std::size_t size) throw()
    {
        return malloc(size);
    }

    void *DoReallocate(void *p, std::size_t size) throw()
    {
        return realloc(p, size);
    }

    void DoDeallocate(void *p, std::size_t size) throw()
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

// override new/delete operators using our allocators
void *operator new(size_t size, StdAllocator &alloc)
{
    return alloc.Allocate(size);
}

void *operator new[](size_t size, StdAllocator &alloc)
{
    return alloc.Allocate(size);
}

void *operator new(size_t size, StdAllocator &alloc,
                   const std::nothrow_t &)
{
    return alloc.Allocate(size, std::nothrow);
}

void *operator new[](size_t size, StdAllocator &alloc,
                     const std::nothrow_t &)
{
    return alloc.Allocate(size, std::nothrow);
}

void operator delete(void *p, StdAllocator &alloc)
{
    return alloc.Deallocate(p);
}

void operator delete[](void *p, StdAllocator &alloc)
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
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

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
    
    pointer reallocate(size_type n, void *addr = 0)
    {
        return static_cast<pointer>(
            Alloc::GetInstance().Reallocate(addr, n*sizeof(T)));
    }

    void deallocate(pointer p, size_type n)
    {
        Alloc::GetInstance().Deallocate(p);
    }

    size_type max_size() const throw()
    {
        return std::numeric_limits<std::size_t>::max() / sizeof(T);
    }

    void construct(pointer p, const_reference val)
    {
        new ((void *) p) value_type(val);
    }

    void destroy(pointer p)
    {
        p->~value_type();
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

#endif  // ALLOC_H
