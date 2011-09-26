/* -*- C++ -*-
 *
 * mmf.h -- Matrix Market Format routines
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011,      Vasileios Karakasis
 * Copyright (C) 2011,      Theodorοs Gκountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef MMF_H__
#define MMF_H__

#include "spm.h"

#include <iostream>
#include <iterator>

namespace csx {

void getMmfHeader(const char *mmf_file,
                  uint64_t &nrows, uint64_t &ncols, uint64_t &nnz);
void getMmfHeader(std::istream &in,
                  uint64_t &nrows, uint64_t &ncols, uint64_t &nnz);


class MMF
{
public:
    uint64_t nrows, ncols, nnz;

    // initialization
    MMF(std::istream &in);
    
    // get next element (false if end)
    bool next(uint64_t &y, uint64_t &x, double &val);
    
    // CooElem iterator
    class iterator;
    iterator begin();
    iterator end();

private:
    std::istream &in_;
};

class MMF::iterator : public std::iterator<std::forward_iterator_tag, CooElem>
{
public:
    iterator() {}
    
    iterator(MMF *mmf, uint64_t cnt) : mmf_(mmf), cnt_(cnt)
    {
        // this is the initializer
        if (cnt_ == 0) {
            this->DoSet();
        }
    }

    bool operator==(const iterator &i)
    {
        //std::cout << "me: " << mmf_ << " " << cnt_
        //          << " i: " << i.mmf << " " << i.cnt << "\n";
        return (mmf_ == i.mmf_) && (cnt_ == i.cnt_);
    }

    bool operator!=(const iterator &i)
    {
        return !(*this == i);
    }

    void operator++()
    {
        ++cnt_;
        this->DoSet();
    }

    CooElem operator*()
    {
        if (!valid_) {
            std::cout << "Requesting dereference, but mmf ended\n"
                      << "cnt: " << cnt_ << std::endl;
            assert(false);
        }
        
        assert(valid_);
        return elem_;
    }

private:
    void DoSet()
    {
        valid_ = mmf_->next(elem_.y, elem_.x, elem_.val);
    }

    MMF *mmf_;
    uint64_t cnt_;
    CooElem elem_;
    bool valid_;
};

MMF::iterator MMF::begin()
{
    return MMF::iterator(this, 0);
}

MMF::iterator MMF::end()
{
    return MMF::iterator(this, this->nnz);
}

} // csx namespace end

#endif  // MMF_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
