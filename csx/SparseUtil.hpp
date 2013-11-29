/* -*- C++ -*-
 *
 * SparseUtil.hpp --  Essentials of SparseInternal/SparsePartition.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSE_UTIL_HPP
#define SPARSE_UTIL_HPP

#include "Encodings.hpp"

#include <ostream>
#include <iomanip>
#include <inttypes.h>
#include <map>
#include <stdint.h>
#include <boost/function.hpp>

using namespace std;

namespace csx {

inline int gcd(int i, int j)
{
    if (j == 0)
        return i;
    return gcd(j, i%j);
}

inline int lcm(int i, int j)
{
    if (i >= j)
        return i * j / gcd(i,j);
    else
        return i * j / gcd(j,i);
}

/**
 *  Coordinate element, i.e., holds row and column information.
 */
template<typename IndexType, typename ValueType>
struct CooElem {
    IndexType row;        ///< the row index
    IndexType col;        ///< the column index
    union {
        ValueType val;   ///< the value of the element
        ValueType *vals; ///< the value of the elements, if Elem refers to
                         ///  an encoded pattern
    };

    CooElem(IndexType r, IndexType c, ValueType v)
        : row(r), col(c), val(v)
    { }

    CooElem(const CooElem &other)
    {
        row = other.row;
        col = other.col;
        val = other.val;
    }

    CooElem()
    {
        row = col = 0;
        val = 0;
    }

    virtual ~CooElem()
    { }
};

/**
 *  "Less" Sorting Functor for coordinate elements. Can be used by std::sort.
 */ 
template<typename IndexType, typename ValueType>
struct CooElemSorter {
    bool operator() (const CooElem<IndexType, ValueType> &lhs,
                     const CooElem<IndexType, ValueType> &rhs) const
    {
        if (lhs.row < rhs.row)
            return true;

        if (lhs.row > rhs.row)
            return false;

        if (lhs.col < rhs.col)
            return true;

        if (lhs.col > rhs.col)
            return false;

        return false;
    }
};

template<typename IndexType, typename ValueType>
struct VerticalSorter {
    public:
        bool operator() (const CooElem<IndexType, ValueType> &lhs,
                         const CooElem<IndexType, ValueType> &rhs) const
        {
            if (lhs.col < rhs.col)
                return true;

            if (lhs.col > rhs.col)
                return false;

            if (lhs.row < rhs.row)
                return true;

            if (lhs.row > rhs.row)
                return false;

            return false;
        }
};

template<typename IndexType, typename ValueType>
struct DiagonalSorter {
    public:
        bool operator() (const CooElem<IndexType, ValueType> &lhs,
                         const CooElem<IndexType, ValueType> &rhs) const
        {
            if ((lhs.col - lhs.row) < (rhs.col - rhs.row)) return true;
            if ((lhs.col - lhs.row) > (rhs.col - rhs.row)) return false;

            // if (min(lhs.row, lhs.col) < min(rhs.row, rhs.col)) return true;
            // if (min(lhs.row, lhs.col) > min(rhs.row, rhs.col)) return false;

            return false;
        }
};

template<typename IndexType, typename ValueType>
struct RevDiagonalSorter {
    public:
        RevDiagonalSorter(int ncols) : ncols(ncols) {}
        bool operator() (const CooElem<IndexType, ValueType> &lhs,
                         const CooElem<IndexType, ValueType> &rhs) const
        {
            if ((lhs.col + lhs.row) < (rhs.col + rhs.row)) return true;
            if ((lhs.col + lhs.row) > (rhs.col + rhs.row)) return false;

            IndexType tmp_col_lhs = ((lhs.col + lhs.row - 1) <= ncols)
                ? lhs.row : ncols + 1 - lhs.col;
            IndexType tmp_col_rhs = ((rhs.col + rhs.row - 1) <= ncols)
                ? rhs.row : ncols + 1 - rhs.col;
            if (tmp_col_lhs < tmp_col_rhs) return true;
            if (tmp_col_lhs > tmp_col_rhs) return false;

            return false;
        }
    private:
        int ncols;
};

template<typename IndexType, typename ValueType>
struct BlockRowSorter {
    public:
        BlockRowSorter(int r) : r(r) {}
        bool operator() (const CooElem<IndexType, ValueType> &lhs,
                         const CooElem<IndexType, ValueType> &rhs) const
        {
            if (((lhs.row - 1) / r) < ((rhs.row - 1) / r)) return true;
            if (((lhs.row - 1) / r) > ((rhs.row - 1) / r)) return false;
                
            if (((lhs.row - 1) % r + r * (lhs.col - 1)) <
                ((rhs.row - 1) % r + r * (rhs.col - 1))) return true;
            if (((lhs.row - 1) % r + r * (lhs.col - 1)) >
                ((rhs.row - 1) % r + r * (rhs.col - 1))) return false;

            return false;
        }
    private:
        int r;
};

template<typename IndexType, typename ValueType>
struct BlockColSorter {
    public:
        BlockColSorter(int c) : c(c) {}
        bool operator() (const CooElem<IndexType, ValueType> &lhs,
                         const CooElem<IndexType, ValueType> &rhs) const
        {
            if (((lhs.col - 1) / c) < ((rhs.col - 1) / c)) return true;
            if (((lhs.col - 1) / c) > ((rhs.col - 1) / c)) return false;
                
            if (((lhs.col - 1) % c + c * (lhs.row - 1)) <
                ((rhs.col - 1) % c + c * (rhs.row - 1))) return true;
            if (((lhs.col - 1) % c + c * (lhs.row - 1)) >
                ((rhs.col - 1) % c + c * (rhs.row - 1))) return false;

            return false;
        }  
    private:
        int c;
};


/**
 * Compares two coordinate elements. This function imposes a lexicographical
 * order in the elements of the matrix.
 *
 * @param p0 the first CooElem to compare
 * @param p1 the second CooElem to compare
 * @return   1 if <tt>p0</tt> succeeds <tt>p1</tt> in lexicographical order,
 *          -1 if <tt>p0</tt> preceeds <tt>p1</tt> in lexicographical order,
 *           0 if 'p0' and 'p1' refer to the same matrix element.
 *  @see CooElem
 */
template<typename IndexType, typename ValueType>
inline int CooCmp(const CooElem<IndexType, ValueType> &p0,
                  const CooElem<IndexType, ValueType> &p1)
{
    int32_t ret;// here is the problem

    ret = p0.row - p1.row;
    if (ret == 0)
        ret = p0.col - p1.col;

    if (ret > 0)
        return 1;
    else if (ret < 0)
        return -1;
    else
        return 0;
}

/**
 *  The CSX pattern id, which will be encoded in uflags.
 *
 *  @see ctl_ll.h
 */
#define CSX_PID_OFFSET 10000L
class DeltaRLE
{
public:
    DeltaRLE(size_t size, size_t delta, Encoding::Type type)
        : size_(size), delta_(delta)
    {
        type_ = type;
    }

    virtual DeltaRLE *Clone() const
    {
        return new DeltaRLE(*this);
    }

    virtual ~DeltaRLE() {};

    /**
     *  Calculates the size of patterns in bytes.
     *
     *  @return size of the selected object of class DeltaRLE.
     */
    virtual size_t GetSize() const
    {
        return size_;
    }

    /**
     *  Gets the type of this RLE pattern.
     *
     *  @return the type of this RLE pattern.
     */
    virtual Encoding::Type GetType() const
    {
        return type_;
    }

    virtual size_t GetDelta() const
    {
        return delta_;
    }

    virtual size_t GetOtherDim() const
    {
        return 0;
    }
    /**
     *  Retrieves the id of this pattern.
     *
     *  @return the id of this pattern.
     */
    virtual size_t GetPatternId() const
    {
        assert(type_ > Encoding::None && type_ < Encoding::BlockRowMin &&
               "not a drle type");
        return type_ * CSX_PID_OFFSET + delta_;
    }

    /**
     *  Calculates the difference between the current and the next element.
     *
     *  @param  order a type of pattern.
     *  @return       number of columns to procceed in order to go to the next
     *                CSX element.
     */
    virtual size_t ColIncrease(Encoding::Type order) const
    {
        size_t ret = (order == type_) ? (size_ * delta_) : 1;
        return ret;
    }

    /**
     *  Calculates the difference between the current and the next element.
     *
     *  @param order a type of pattern.
     *  @param jmp   a number of columns to jump.
     *  @return      number of columns to procceed in order to go to the next
     *               CSX element.
     */
    virtual size_t ColIncreaseJmp(Encoding::Type order, size_t jmp) const
    {
        size_t ret = jmp;
        if (order == type_)
            ret += ((size_ - 1) * delta_);

        return ret;
    }

    /**
     *  Prints the data of this pattern on the specified output stream.
     *
     *  @param out the output stream.
     */
    virtual std::ostream &PrintOn(std::ostream &out) const
    {
        out << "drle: size=" << size_ << " len=" << delta_
            << " type=" << type_;
        return out;
    }

    /**
     *  Handles the elements within pattern. Useful for debugging.
     */
    class Generator
    {
    public:
        Generator(CooElem<size_t, double> start, DeltaRLE *rle)
            : start_(start), rle_(rle), nr_(0) { }

        /**
         *  Checks if the handler is in the last element of the pattern.
         *
         *  @return true if the handler is in the last element and false
         *          otherwise
         */
        virtual bool IsEmpty() const {
            return (nr_ == rle_->size_);
        }

        /**
         *  Finds the next element within pattern.
         *
         *  @return next element
         */
        virtual CooElem<size_t, double> Next() {
            CooElem<size_t, double> ret(start_);
            assert(nr_ <= rle_->size_ && "out of pattern");
            ret.col += nr_ * rle_->delta_;
            nr_ += 1;
            return ret;
        }

    private:
        CooElem<size_t, double> start_;
        DeltaRLE *rle_;
        size_t nr_;
    };

    /**
     *  Creates an object of Generator class.
     *
     *  @param  start starting element
     *  @return an object of Generator class.
     */
    virtual Generator *generator(CooElem<size_t, double> start)
    {
        Generator *g = new DeltaRLE::Generator(start, this);
        return g;
    }

    /**
     *  Finds the column of the next non-zero element within pattern.
     *
     *  @return column of the next element of pattern.
     */
    virtual size_t GetNextCol(size_t x0) const
    {
        return x0 + delta_;
    }

    /**
     *  Stats for a specific type of pattern of sparse matrix. More
     *  specifically, it just holds the number of non-zero elements and
     *  patterns that adhere to this specific type of pattern.
     */
    struct StatsVal
    {
        size_t nnz;
        size_t npatterns;

        /**
         *  Update the stats for every new instance of this type of pattern.
         *
         *  @param new_vals stats of an instance of a pattern.
         */
        void Update(const StatsVal& new_vals)
        {
            this->nnz += new_vals.nnz;
            this->npatterns += new_vals.npatterns;
        }
    };

    // Stats[0] -> delta values (encoded with this DRLE)
    typedef std::map<size_t, DeltaRLE::StatsVal> Stats;

protected:
    size_t size_, delta_;

private:
    Encoding::Type type_;
};

/**
 *  A block RLE pattern. It is essentially a DeltaRLE pattern with
 *  <tt>delta=1</tt> augmented with one more field for the second dimension
 *  of the block pattern.
 */
class BlockRLE : public DeltaRLE {
public:
    BlockRLE(size_t size, size_t other_dim, Encoding &enc)
        : DeltaRLE(size, 1, (assert(enc.GetBlockAlignment()), enc.GetType()))
    {
        other_dim_ = other_dim;
    }

    /**
     *  Finds the second dimension of blocks (first is revealed from type).
     *
     *  @return the size of second dimension.
     */
    virtual size_t GetOtherDim() const
    {
        return other_dim_;
    }

    virtual unsigned long GetPatternId() const
    {
        return CSX_PID_OFFSET * GetType() + other_dim_;
    }

    virtual BlockRLE *Clone() const
    {
        return new BlockRLE(*this);//memory leak
    }

private:
    size_t other_dim_;
};

/**
 *  A generic sparse matrix pattern.
 */
class SpmPattern
{
public:
    DeltaRLE *pattern;
    bool in_pattern;    // only for statistics
    bool pattern_start; // only for statistics (pattern start)

    SpmPattern() : pattern(NULL), in_pattern(false), pattern_start(false) { }
    SpmPattern(const SpmPattern &spm_p)
    {
        DeltaRLE *p;
        this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->Clone();
        this->in_pattern = spm_p.in_pattern;
        this->pattern_start = spm_p.pattern_start;
    }

    virtual ~SpmPattern()
    {
        if (this->pattern != NULL)
            delete this->pattern;
    }

    SpmPattern& operator=(const SpmPattern &spm_p)
    {
        if (this == &spm_p)
            return *this;

        if (this->pattern != NULL)
            delete this->pattern;

        DeltaRLE *p;

        this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->Clone();
        this->in_pattern = spm_p.in_pattern;
        this->pattern_start = spm_p.pattern_start;
        return *this;
    }
};

/**
 *  A generic sparse matrix coordinate element that can also be a pattern.
 */
template<typename IndexType, typename ValueType> 
struct Elem : public CooElem<IndexType, ValueType>, public SpmPattern
{
    Elem(const CooElem<IndexType, ValueType> &other)
        : CooElem<IndexType, ValueType>(other)
    {
        pattern = 0;
    }

    Elem(const Elem<IndexType, ValueType> &other)
        : CooElem<IndexType, ValueType>(other.row, other.col, other.val)
    {
        in_pattern = other.in_pattern;
        pattern_start = other.pattern_start;
        pattern = (other.pattern) ? (other.pattern)->Clone() : NULL;
    }

    Elem()
        : CooElem<IndexType, ValueType>(0, 0, 0)
    {
        pattern = 0;
    }
};


template<typename IndexType, typename ValueType>
inline bool elem_cmp_less(const Elem<IndexType, ValueType> &e0,
                          const Elem<IndexType, ValueType> &e1)
{
    int ret;
    ret = CooCmp(static_cast<CooElem<IndexType, ValueType> >(e0),
                 static_cast<CooElem<IndexType, ValueType> >(e1));
    return (ret < 0);
}

/**
 *  Fills a sparse matrix element (<tt>dst</tt>) from a source element
 *  (<tt>src</tt>) element.
 *
 *  @param src  element to copy from.
 *  @param dst  element to fill.
 */
template<typename IndexType, typename ValueType>
void MakeRowElem(const CooElem<IndexType, ValueType> &p,
                 Elem<IndexType, ValueType> *ret)
{
    ret->row = p.row;
    ret->col = p.col;
    ret->val = p.val;
    ret->pattern = NULL;
}

template<typename IndexType, typename ValueType>
void MakeRowElem(const Elem<IndexType, ValueType> &p,
                 Elem<IndexType, ValueType> *ret)
{
    ret->row = p.row;
    ret->col = p.col;
    ret->val = p.val;
    ret->in_pattern = p.in_pattern;
    ret->pattern_start = p.pattern_start;
    ret->pattern = (p.pattern == NULL) ? NULL : (p.pattern)->Clone();
}

/**
 * Stream operations
 */
std::ostream &operator<<(std::ostream &os, const DeltaRLE::StatsVal &stats)
{
    os << "nnz: " << stats.nnz;
    return os;
}

std::ostream &operator<<(std::ostream &os, const DeltaRLE &p)
{
    os << " (";
    p.PrintOn(os);
    os << " type:" << p.GetType() << ") ";
    return os;
}

template<typename IndexType, typename ValueType>
std::ostream &operator<<(std::ostream &out, CooElem<IndexType, ValueType> p)
{
    out << "(" << p.row << "," << p.col << ")";
    return out;
}

template<typename IndexType, typename ValueType>
std::ostream &operator<<(std::ostream &out, const Elem<IndexType, ValueType> e)
{
    out << static_cast<CooElem<IndexType, ValueType> >(e);
    if (e.pattern != NULL) {
        out << "->[" << *(e.pattern) << "]";
        out << " vals:{ ";
        for (size_t i = 0; i < e.pattern->GetSize(); i++)
            out << e.vals[i] << " ";
            
        out << "}";
    } else {
        out << "v:" << e.val;
    }

    return out;
}

// template<typename IndexType, typename ValueType>
// std::ostream &operator<<(std::ostream &out,
//                          const Elem<IndexType, ValueType> &elem)
// {
//     out << "row:" << elem.row;
//     out << "col:" << elem.col;
//     out << " in_pattern:" << elem.in_pattern;
//     out << " pattern_start:" << elem.pattern_start;
//     if (elem.pattern) {
//         out << *(elem.pattern);
//         out << " vals:{ ";
//         for (int i = 0; i < elem.pattern->GetSize(); i++)
//             out << elem.vals[i] << " ";
            
//         out << "}";
//     } else {
//         out << " v:" << elem.val;
//     }

//     return out;
// }

template<typename IndexType, typename ValueType>
struct TransformFnType {
    typedef boost::function<void (CooElem<IndexType, ValueType> &p)> TransformFn;
};

}   // namespace csx

#endif  // SPARSE_UTIL_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
