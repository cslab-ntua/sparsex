/* -*- C++ -*-
 *
 * spm.h --  Internal representation of sparse matrices.
 *
 * Copyright (C) 2010, Computing Systems Laboratory (CSLab), NTUA.
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SPM_H__
#define CSX_SPM_H__

#include <inttypes.h>
#include <vector>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdint.h>

#include <boost/function.hpp>
#include <boost/foreach.hpp>

#include "dynarray.h"

namespace csx
{

class MMF;

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
} SpmIterOrder;

const SpmIterOrder SpmTypes[] =
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

/**
 * Determines if the input is a block type or not.
 *
 * @param t a type of pattern.
 * @return  block alignment if <tt>t</tt> is a block type, 0 otherwise.
 */
static inline int IsBlockType(SpmIterOrder t)
{
    if (t > BLOCK_TYPE_START && t < BLOCK_COL_START)
        return t - BLOCK_TYPE_START;
    else if (t > BLOCK_COL_START && t < BLOCK_TYPE_END)
        return t - BLOCK_COL_START;
    else
        return 0;
}

/**
 *  Holds column and value or values of CSX elements.
 */
struct RowElem {
    uint64_t x;     ///< the column index
    union {
        double val;     ///< the value of the element
        double *vals;   ///< the value of the elements, if RowElem refers to
                        ///  an encoded pattern
    };
};

/**
 *  Coordinate element, i.e., holds row and column information.
 */
struct CooElem : public RowElem {
    uint64_t y; ///< the row of the element
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
static inline int CooCmp(const CooElem &p0, const CooElem &p1)
{
    int64_t ret;

    ret = p0.y - p1.y;
    if (ret == 0)
        ret = p0.x - p1.x;

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
    DeltaRLE(uint32_t size, uint32_t delta, SpmIterOrder type)
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
    virtual SpmIterOrder GetType() const
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
    virtual long ColIncrease(SpmIterOrder order) const
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
    virtual uint64_t ColIncreaseJmp(SpmIterOrder order, uint64_t jmp) const
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
        Generator(CooElem start, DeltaRLE *rle)
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
        virtual CooElem Next() {
            CooElem ret(start_);
            assert(nr_ <= rle_->size_ && "out of pattern");
            ret.x += nr_ * rle_->delta_;
            nr_ += 1;
            return ret;
        }

    private:
        CooElem start_;
        DeltaRLE *rle_;
        uint64_t nr_;
    };

    /**
     *  Creates an object of Generator class.
     *
     *  @param  start starting element
     *  @return an object of Generator class.
     */
    virtual Generator *generator(CooElem start)
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

    typedef std::map<uint64_t, DeltaRLE::StatsVal> Stats;

protected:
    uint32_t size_, delta_;

private:
    SpmIterOrder type_;
};

/**
 *  A block RLE pattern. It is essentially a DeltaRLE pattern with
 *  <tt>delta=1</tt> augmented with one more field for the second dimension
 *  of the block pattern.
 */
class BlockRLE : public DeltaRLE {
public:
    BlockRLE(uint32_t size, uint32_t other_dim, SpmIterOrder type)
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

    SpmPattern(void) : pattern(NULL) { }
    SpmPattern(const SpmPattern &spm_p)
    {
        DeltaRLE *p;

        this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->Clone();
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
        return *this;
    }
};

/**
 *  A generic sparse matrix coordinate element that can also be a pattern.
 */
class SpmCooElem: public CooElem, public SpmPattern {};

/**
 *  A generic sparse matrix row element that can also be a pattern.
 */
class SpmRowElem: public RowElem, public SpmPattern {};

/**
 *  Fills a sparse matrix row element (<tt>dst</tt>) from a source element
 *  (<tt>src</tt>) element.
 *
 *  @param src  element to copy from.
 *  @param dst  element to fill.
 */
void MakeRowElem(const CooElem &src, SpmRowElem *dst);
void MakeRowElem(const SpmCooElem &src, SpmRowElem *dst);
void MakeRowElem(const SpmRowElem &src, SpmRowElem *dst);

typedef std::iterator<std::forward_iterator_tag, CooElem> SpmPointIter;
typedef boost::function<void (CooElem &p)> TransformFn;

/**
 *  Internal representation of a sparse matrix. This class can also represent a
 *  sub-matrix of a larger one, starting at a specific row of the original
 *  matrix.
 */
class SPM
{
private:
    uint64_t nr_rows_, nr_cols_, nr_nzeros_;
    SpmIterOrder type_;
    SpmRowElem *elems_;
    uint64_t elems_size_;
    uint64_t *rowptr_;
    uint64_t rowptr_size_;
    uint64_t row_start_;    ///< Row of the original matrix, where this
                            ///  sub-matrix starts.
    bool elems_mapped_;

    // These are mainly SPM construction classes, so make them friends
    friend class Builder;
    friend class VirtualBuilder;
    friend class CsxManager;
    friend class DRLE_Manager;

public:
    SPM() : type_(NONE), elems_(NULL), rowptr_(NULL) {}

    ~SPM()
    {
        if (!elems_mapped_ && elems_)
            free(elems_);

        if (rowptr_)
            free(rowptr_);
    };

    uint64_t GetNrRows()
    {
        return rowptr_size_ - 1;
    }

    uint64_t GetNrCols()
    {
        return nr_cols_;
    }

    uint64_t GetNrNonzeros()
    {
        return nr_nzeros_;
    }

    SpmIterOrder GetType()
    {
        return type_;
    }

    SpmRowElem *RowBegin(uint64_t ridx = 0);
    SpmRowElem *RowEnd(uint64_t ridx = 0);

    class Builder;
    class VirtualBuilder;

    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start point iterators start.
     *  @param pnts_end   point iterators end.
     *  @param first_row  first row of the matrix.
     *  @param limit      limit for the number of elements.
     *  @param nr_elems   size for elems allocation.
     *  @param nr_rows    size for rows allocation.
     *  @param SpmBld     point to a builder responsible for allocations.
     *  @return           number of elements allocated.
     */
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit = 0,
                      uint64_t nr_elems = 0, uint64_t nr_rows = 0);
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit,
                      uint64_t nr_elems, uint64_t nr_rows,
                      SPM::Builder *SpmBld);

    /**
     *  Loads a matrix from file.
     *
     *  @param filename name of the file that matrix is kept.
     *  @param in       buffer from which the matrix is taken.
     *  @return         spm class object with the characteristics of the matrix.
     */
    static SPM *LoadMMF(const char *filename);
    static SPM *LoadMMF(std::istream &in=std::cin);

    /**
     *  Loads matrix from a file specifying the number of threads.
     *
     *  @param filename name of the file that matrix is kept.
     *  @param in       buffer from which the matrix is taken.
     *  @param mmf      handler of MMF class.
     *  @param nr       number of threads to be used.
     *  @return         spm class object with the characteristics of the matrix.
     */
    static SPM *LoadMMF_mt(const char *mmf_file, const long nr);
    static SPM *LoadMMF_mt(std::istream &in, const long nr);
    static SPM *LoadMMF_mt(MMF &mmf, const long nr);

    /**
     *  Printing methods.
     */
    void Print(std::ostream &out = std::cout);
    void PrintElems(std::ostream &out = std::cout);
    void PrintRows(std::ostream &out = std::cout);
    void PrintStats(std::ostream &out = std::cout);

    class PntIter;

    PntIter PointsBegin(uint64_t ridx = 0);
    PntIter PointsEnd(uint64_t ridx = 0);

    /**
     *  Transforms matrix to specific type.
     *
     *  @param type type of pattern to transform to.
     *  @param rs   line to start the transformation. If it is zero, it is the
     *              the beggining of matrix.
     *  @param re   line to end the transformation. If it is zero, it is the
     *              the end of matrix.
     */
    void Transform(SpmIterOrder type, uint64_t rs = 0, uint64_t re = 0);

    /**
     *  Finds the reverse transform function. Transforms from <tt>type</tt> to
     *  <tt>HORIZONTAL</tt>.
     *
     *  @param type type of pattern to transform from.
     *  @return     reverse transform function.
     */
    TransformFn GetRevXformFn(SpmIterOrder type);

    /**
     *  Finds the transform function. Transforms from HORIZONTAL to
     *  <tt>type</tt>.
     *
     *  @param type type of pattern to transform to.
     *  @return     transform function.
     */
    TransformFn GetXformFn(SpmIterOrder type);

    /**
     *  Finds the complete transform function. Tranforms from <tt>from</tt> to
     *  <tt>to</tt>.
     *
     *  @param from type of pattern to transform from.
     *  @param to   type of pattern to transform to.
     *  @return     complete transform function.
     */
    TransformFn GetTransformFn(SpmIterOrder from, SpmIterOrder to);

    /**
     *  Extracts a window (sub-matrix) from this matrix. This methods copies the
     *  nonzero elements of the original matrix to the window.
     *
     *  @param rs     the row in this matrix, from where the window will be
     *                extracted.
     *  @param length number of rows to include in the window.
     *  @return       the selected window.
     *  @see GetWindow()
     */
    SPM *ExtractWindow(uint64_t rs, uint64_t length);

    /**
     *  Extracts a window (sub-matrix) from this matrix using a backing
     *  store. This method uses the nonzero elements of this matrix without
     *  performing any copy.
     *
     *  @param rs     the row in this matrix, from where the window will be
     *                extracted.
     *  @param length number of rows to include in the window.
     *  @return       the selected window.
     *  @see ExtractWindow()
     */
    SPM *GetWindow(uint64_t rs, uint64_t length);

    /**
     *  Puts a (possibly processed) window of the matrix back to this matrix.
     *
     *  @param window the (possibly processed) window.
     */
    void PutWindow(const SPM *window);
};

void TestMMF(SPM *spm, const char *mmf_file);

/**
 *  Utility class for building an SPM.
 */
class SPM::Builder
{
public:
    /**
     *  Constructs a Builder object.
     *
     *  @param spm      the matrix to be built.
     *  @param nr_elems an initial guess of the nonzero elements.
     *  @param nr_rows  an initial guess of the number of rows.
     */
    Builder(SPM *spm, uint64_t nr_elems = 0, uint64_t nr_rows = 0);
    virtual ~Builder();

    /**
     *  Allocates a new element in the new matrix.
     *
     *  @return pointer to element allocated.
     */
    virtual SpmRowElem *AllocElem();

    /**
     *  Allocates a number of elements in the new matrix.
     *
     *  @param  nr number of elements to be allocated.
     *  @return    pointer to the first allocated element.
     */
    virtual SpmRowElem *AllocElems(uint64_t nr);

    /**
     *  Counts elements already allocated in the matrix.
     *
     *  @return number of elements allocated.
     */
    virtual uint64_t GetElemsCnt();

    /**
     *  Allocate elements for the next row.
     *
     *  @param rdiff number of empty rows between the current and the new row.
     */
    virtual void NewRow(uint64_t rdiff = 1);

    /**
     *  Finalizes the construction of the new matrix.
     */
    virtual void Finalize();

private:
    SPM *spm_;
    dynarray_t *da_elems_;
    dynarray_t *da_rowptr_;
};

/**
 *  Iterators for moving between elements inside the matrix.
 */
class SPM::PntIter : public std::iterator<std::forward_iterator_tag, CooElem>
{
public:
    PntIter(): spm_(NULL), row_idx_(0), elem_idx_(0) { }
    PntIter(SPM *spm, uint64_t row_idx);

    bool operator==(const PntIter &pi);
    bool operator!=(const PntIter &pi);
    void operator++();
    SpmCooElem operator*();
    std::ostream &operator<<(std::ostream &out);

private:
    SPM *spm_;
    uint64_t row_idx_;
    uint64_t elem_idx_;
};

/**
 * Forward declarations for stream operations
 */
std::ostream &operator<<(std::ostream &os, const DeltaRLE::StatsVal &stats);
std::ostream &operator<<(std::ostream &os, const DeltaRLE &p);
std::ostream &operator<<(std::ostream &out, CooElem p);
std::ostream &operator<<(std::ostream &out, const SpmCooElem e);
std::ostream &operator<<(std::ostream &out, const SpmRowElem &elem);
}   // end of csx namespace

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
