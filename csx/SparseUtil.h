/* -*- C++ -*-
 *
 * SparseUtil.h --  Essentials of SparseInternal/SparsePartition.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSEUTIL_H__
#define SPARSEUTIL_H__

#include <ostream>
#include <iomanip>
#include <inttypes.h>
#include <map>
#include <stdint.h>
#include <boost/function.hpp>

namespace csx {

#define STRINGIFY__(s) #s
#define STRINGIFY(s)  STRINGIFY__(s)
#define BLOCK_ROW_TYPE_NAME(r)  BLOCK_R ## r
#define BLOCK_COL_TYPE_NAME(c)  BLOCK_C ## c

///< Supported encoding patterns for CSX.
typedef enum
{
    NONE = 0,
    HORIZONTAL,
    VERTICAL,
    DIAGONAL,
    REV_DIAGONAL,
    BLOCK_TYPE_START,
    BLOCK_ROW_TYPE_NAME(1),
    BLOCK_ROW_TYPE_NAME(2),
    BLOCK_ROW_TYPE_NAME(3),
    BLOCK_ROW_TYPE_NAME(4),
    BLOCK_ROW_TYPE_NAME(5),
    BLOCK_ROW_TYPE_NAME(6),
    BLOCK_ROW_TYPE_NAME(7),
    BLOCK_ROW_TYPE_NAME(8),
    BLOCK_COL_START,
    BLOCK_COL_TYPE_NAME(1),
    BLOCK_COL_TYPE_NAME(2),
    BLOCK_COL_TYPE_NAME(3),
    BLOCK_COL_TYPE_NAME(4),
    BLOCK_COL_TYPE_NAME(5),
    BLOCK_COL_TYPE_NAME(6),
    BLOCK_COL_TYPE_NAME(7),
    BLOCK_COL_TYPE_NAME(8),
    BLOCK_TYPE_END,
    XFORM_MAX
} IterOrder;

const IterOrder SpmTypes[] =
{
    NONE,
    HORIZONTAL,
    VERTICAL,
    DIAGONAL,
    REV_DIAGONAL,
    BLOCK_TYPE_START,
    BLOCK_ROW_TYPE_NAME(1),
    BLOCK_ROW_TYPE_NAME(2),
    BLOCK_ROW_TYPE_NAME(3),
    BLOCK_ROW_TYPE_NAME(4),
    BLOCK_ROW_TYPE_NAME(5),
    BLOCK_ROW_TYPE_NAME(6),
    BLOCK_ROW_TYPE_NAME(7),
    BLOCK_ROW_TYPE_NAME(8),
    BLOCK_COL_START,
    BLOCK_COL_TYPE_NAME(1),
    BLOCK_COL_TYPE_NAME(2),
    BLOCK_COL_TYPE_NAME(3),
    BLOCK_COL_TYPE_NAME(4),
    BLOCK_COL_TYPE_NAME(5),
    BLOCK_COL_TYPE_NAME(6),
    BLOCK_COL_TYPE_NAME(7),
    BLOCK_COL_TYPE_NAME(8),
    BLOCK_TYPE_END,
    XFORM_MAX
};

///< Names of the encoding types
const char *SpmTypesNames[] =
{
    "__NONE__",
    "HORIZONTAL",
    "VERTICAL",
    "DIAGONAL",
    "REV_DIAGONAL",
    "__BLOCK_TYPE_START__",
    STRINGIFY(BLOCK_ROW_TYPE_NAME(1)),
    STRINGIFY(BLOCK_ROW_TYPE_NAME(2)),
    STRINGIFY(BLOCK_ROW_TYPE_NAME(3)),
    STRINGIFY(BLOCK_ROW_TYPE_NAME(4)),
    STRINGIFY(BLOCK_ROW_TYPE_NAME(5)),
    STRINGIFY(BLOCK_ROW_TYPE_NAME(6)),
    STRINGIFY(BLOCK_ROW_TYPE_NAME(7)),
    STRINGIFY(BLOCK_ROW_TYPE_NAME(8)),
    "__BLOCK_COL_START__",
    STRINGIFY(BLOCK_COL_TYPE_NAME(1)),
    STRINGIFY(BLOCK_COL_TYPE_NAME(2)),
    STRINGIFY(BLOCK_COL_TYPE_NAME(3)),
    STRINGIFY(BLOCK_COL_TYPE_NAME(4)),
    STRINGIFY(BLOCK_COL_TYPE_NAME(5)),
    STRINGIFY(BLOCK_COL_TYPE_NAME(6)),
    STRINGIFY(BLOCK_COL_TYPE_NAME(7)),
    STRINGIFY(BLOCK_COL_TYPE_NAME(8)),
    "__BLOCK_TYPE_END__",
    "__XFORM_MAX__"
};

inline bool is_row_block(IterOrder t)
{
    if (t > BLOCK_TYPE_START && t < BLOCK_COL_START)
        return true;
    return false;
}

inline bool is_col_block(IterOrder t)
{
    if (t > BLOCK_COL_START && t < BLOCK_TYPE_END)
        return true;
    return false;
}

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
 * Determines if the input is a block type or not.
 *
 * @param t a type of pattern.
 * @return  block alignment if <tt>t</tt> is a block type, 0 otherwise.
 */
inline int IsBlockType(IterOrder t)
{
    if (t > BLOCK_TYPE_START && t < BLOCK_COL_START)
        return t - BLOCK_TYPE_START;
    else if (t > BLOCK_COL_START && t < BLOCK_TYPE_END)
        return t - BLOCK_COL_START;
    else
        return 0;
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
};

/**
 *  "Less" Sorting Functor for coordinate elements. Can be used by std::sort.
 */ 
template<typename IndexType, typename ValueType>
struct CooElemSorter {
    public:
        bool operator() (const CooElem<IndexType, ValueType> &lhs,
                         const CooElem<IndexType, ValueType> &rhs) const
        {
            if (lhs.row < rhs.row) return true;
            if (lhs.row > rhs.row) return false;

            if (lhs.col < rhs.col) return true;
            if (lhs.col > rhs.col) return false;

            return false;
        }
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
    int64_t ret;

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
    DeltaRLE(uint32_t size, uint32_t delta, IterOrder type)
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
    virtual long GetSize() const
    {
        return size_;
    }

    /**
     *  Gets the type of this RLE pattern.
     *
     *  @return the type of this RLE pattern.
     */
    virtual IterOrder GetType() const
    {
        return type_;
    }

    /**
     *  Retrieves the id of this pattern.
     *
     *  @return the id of this pattern.
     */
    virtual long GetPatternId() const
    {
        assert(type_ > NONE && type_ < BLOCK_TYPE_START &&
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
    virtual long ColIncrease(IterOrder order) const
    {
        long ret = (order == type_) ? (size_ * delta_) : 1;
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
    virtual uint64_t ColIncreaseJmp(IterOrder order, uint64_t jmp) const
    {
        long ret = jmp;
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
        Generator(CooElem<uint64_t, double> start, DeltaRLE *rle)
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
        virtual CooElem<uint64_t, double> Next() {
            CooElem<uint64_t, double> ret(start_);
            assert(nr_ <= rle_->size_ && "out of pattern");
            ret.col += nr_ * rle_->delta_;
            nr_ += 1;
            return ret;
        }

    private:
        CooElem<uint64_t, double> start_;
        DeltaRLE *rle_;
        uint64_t nr_;
    };

    /**
     *  Creates an object of Generator class.
     *
     *  @param  start starting element
     *  @return an object of Generator class.
     */
    virtual Generator *generator(CooElem<uint64_t, double> start)
    {
        Generator *g = new DeltaRLE::Generator(start, this);
        return g;
    }

    /**
     *  Finds the column of the next non-zero element within pattern.
     *
     *  @return column of the next element of pattern.
     */
    virtual uint64_t GetNextCol(uint64_t x0) const
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
        uint64_t nnz;
        long npatterns;

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
    typedef std::map<uint64_t, DeltaRLE::StatsVal> Stats;

protected:
    uint32_t size_, delta_;

private:
    IterOrder type_;
};

/**
 *  A block RLE pattern. It is essentially a DeltaRLE pattern with
 *  <tt>delta=1</tt> augmented with one more field for the second dimension
 *  of the block pattern.
 */
class BlockRLE : public DeltaRLE {
public:
    BlockRLE(uint32_t size, uint32_t other_dim, IterOrder type)
        : DeltaRLE(size, 1, (assert(IsBlockType(type)), type))
    {
        other_dim_ = other_dim;
    }

    /**
     *  Finds the second dimension of blocks (first is revealed from type).
     *
     *  @return the size of second dimension.
     */
    virtual uint32_t GetOtherDim() const
    {
        return other_dim_;
    }

    virtual long GetPatternId() const
    {
        return CSX_PID_OFFSET * GetType() + other_dim_;
    }

    virtual BlockRLE *Clone() const
    {
        return new BlockRLE(*this);
    }

private:
    uint32_t other_dim_;
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

    SpmPattern(void) : pattern(NULL), in_pattern(false), pattern_start(false) { }
    SpmPattern(const SpmPattern &spm_p)
    {
        DeltaRLE *p;
        this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->Clone();
        this->in_pattern = spm_p.in_pattern;
        this->pattern_start = spm_p.pattern_start;
    }

    ~SpmPattern()
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
struct Elem : public CooElem<IndexType, ValueType>, public SpmPattern {};


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
    out << "(" << std::setw(2) << p.row << "," << std::setw(2) << p.col << ")";
    return out;
}

template<typename IndexType, typename ValueType>
std::ostream &operator<<(std::ostream &out, const Elem<IndexType, ValueType> e)
{
    out << static_cast<CooElem<IndexType, ValueType> >(e);
    if (e.pattern != NULL) {
        out << "->[" << *(e.pattern) << "]";
        out << " vals:{ ";
        for (int i = 0; i < e.pattern->GetSize(); i++)
            out << e.vals[i] << " ";
            
        out << "}\n";
    } else {
        out << "v: " << e.val;
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

#endif  // SPARSEUTIL_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
