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

extern "C" {
    #include "dynarray.h"
}

namespace csx
{

class MMF;

#define STRINGIFY__(s) #s
#define STRINGIFY(s)  STRINGIFY__(s)
#define BLOCK_ROW_TYPE_NAME(r)  BLOCK_R ## r
#define BLOCK_COL_TYPE_NAME(c)  BLOCK_C ## c

/**.
 *  Order patterns used in CSX and match them with numbers.
 */
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
 * @param t a type of pattern.
 * @return  block alignment if 't' is a block type, 0 otherwise.
 */
static inline int IsBlockType(SpmIterOrder t)
{
    if (t > BLOCK_TYPE_START && t < BLOCK_COL_START)
        return t - BLOCK_TYPE_START;
    else if (t > BLOCK_COL_START && t < BLOCK_TYPE_END)
        return t-BLOCK_COL_START;
    else
        return 0;
}

/**
 *  Holds column and value or values of CSX elements.
 */
struct RowElem
{
    uint64_t x;
    union
    {
        double val;
        double *vals;
    };
};

/**
 *  Expansion of RowElem struct which additionally holds row of CSX elements.
 */
struct CooElem : public RowElem
{
    uint64_t y;
};

/**
 * Compare column of two elements.
 * @param p0 an element with row, column and value or values.
 * @param p1 an element with row, column and value or values.
 * @return   1 if 'p0' points to an element after the element that 'p1' points,
 *          -1 if 'p0' points to an element before the element that 'p1' points,
 *           0 if 'p0' and 'p1' point to the same element.
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
 *  Describes pattern of non-zero elements.
 */
#define PID_OFFSET 10000L
class DeltaRLE
{
public:
    SpmIterOrder type;
    uint32_t size, delta;
    
    DeltaRLE(uint32_t _size, uint32_t _delta, SpmIterOrder _type)
            : size(_size), delta(_delta)
    {
        this->type = _type;
    }

    virtual DeltaRLE *Clone() const
    {
        return new DeltaRLE(*this);
    }
    
    virtual ~DeltaRLE() {};

    /**
     *  Calculates the size of patterns in bytes.
     *  @return size of the selected object of class DeltaRLE.
     */
    virtual long GetSize() const
    {
        return this->size;
    }
    
    /**
     *  Generates id for patterns.
     *  @return unique id of the type of the selected object of class DeltaRLE.
     */
    virtual long GetPatId() const
    {
        assert(this->type > NONE && this->type < BLOCK_TYPE_START);
        return type*PID_OFFSET + this->delta;
    }

    /**
     *  Calculates the difference between the current and the next element.
     *  @param  order a type of pattern.
     *  @return       number of columns to procceed in order to go to the next
     *                CSX element.
     */
    virtual long ColIncrease(SpmIterOrder order) const
    {
        long    ret;
        
        ret = (order == this->type) ? (this->size*this->delta) : 1;
        return ret;
    }

    /**
     *  Calculates the difference between the current and the next element.
     *  @param order a type of pattern.
     *  @param jmp   a number of columns to jump.
     *  @return       number of columns to procceed in order to go to the next
     *                CSX element.
     */
    virtual uint64_t ColIncreaseJmp(SpmIterOrder order, uint64_t jmp) const
    {
        long ret;
        
        ret = jmp;
        if (order == this->type)
            ret += ((this->size-1)*this->delta);
        return ret;
    }

    /**
     *  Prints the data of a pattern.
     *  @param out output channel.
     */
    virtual std::ostream &PrintOn(std::ostream &out) const
    {
        out << "drle: size=" << this->size << " len=" << this->delta << " type="
            << this->type;
        return out;
    }

    /**
     *  Handles the elements within pattern. Useful class for debugging.
     */
    class Generator
    {
    public:
    	CooElem start;
	DeltaRLE *rle;
	long nr;

        Generator(CooElem start_, DeltaRLE *rle_): 
                  start(start_), rle(rle_), nr(0) { }

        /**
         *  Checks if the handler is in the last element of the pattern.
         *  @return true if the handler is in the last element and false
         *          otherwise
         */
	virtual bool IsEmpty() const {
		return (this->nr == this->rle->size);
	}
	
	/**
         *  Finds the next element within pattern.
         *  @return next element
         */
	virtual CooElem Next() {
		CooElem ret(start);
		assert(this->nr <= this->rle->size);
		ret.x += (this->nr)*this->rle->delta;
		this->nr += 1;
		return ret;
	}
    };
    
    virtual Generator *generator(CooElem start)
    {
        Generator *g;
        g = new DeltaRLE::Generator(start, this);
        return g;
    }

    /**
     *  Find the column of the next non-zero element within pattern.
     *  @return column of the next element of pattern.
     */
    virtual uint64_t GetNextCol(uint64_t x0) const
    {
        return (x0 + this->delta);
    }

    /**
     *  Stats for a specific type of pattern of sparse matrix. More
     *  specifically, it just holds the number of non-zero elements and
     *  patterns that adhere to this specific type of pattern.
     */
    class StatsVal
    {
    public:
        uint64_t nnz;
        long npatterns;

        /**
         *  Update the stats for every new instance of this type of pattern.
         *  @param new_vals stats of an instance of a pattern.
         */
        void Update(const StatsVal& new_vals)
        {
            this->nnz += new_vals.nnz;
            this->npatterns += new_vals.npatterns;
        }
    };
    
    typedef std::map<uint64_t, DeltaRLE::StatsVal> Stats;
};

/**
 *  Expansion of DeltaRLE to support block patterns.
 */
class BlockRLE : public DeltaRLE {
public:
    uint32_t    other_dim;
    
    BlockRLE(uint32_t size_, uint32_t other_dim_, SpmIterOrder type_)
            : DeltaRLE(size_, 1, (assert(IsBlockType(type_)), type_))
    {
        this->other_dim = other_dim_;
    }
    
    /**
     *  Finds the second dimension of blocks (first is revealed from type).
     *  @return the size of second dimension.
     */
    virtual uint32_t getOtherDim() const
    {
        return this->other_dim;
    }
    
    virtual long GetPatId() const
    {
        return PID_OFFSET*this->type + this->other_dim;
    }
    
    virtual BlockRLE *Clone() const
    {
        return new BlockRLE(*this);
    }
};

/**
 *  Update of DeltaRLE class with two constructors, a destructor and pointer
 *  operator.
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
 *  Expansion of CooElem which holds pattern info about the CSX elements.
 */
class SpmCooElem: public CooElem, public SpmPattern {};

/**
 *  Expansion of RowElem which holds pattern info about the CSX elements.
 */
class SpmRowElem: public RowElem, public SpmPattern {};

/**
 *  Copies elements.
 *  @param p    element of struct CooElem, SpmCooElem or SpmRowElem used to
 *              copy from.
 *  @param ret  element of struct SpmRowElem used to write to.
 */
void MakeRowElem(const CooElem &p, SpmRowElem *ret);
void MakeRowElem(const SpmCooElem &p, SpmRowElem *ret);
void MakeRowElem(const SpmRowElem &p, SpmRowElem *ret);

typedef std::iterator<std::forward_iterator_tag, CooElem> SpmPointIter;
typedef boost::function<void (CooElem &p)> TransformFn;

/**
 *  Sparse Matrix's sub-matrix characteristics. Each thread handles one
 *  sub-matrix (one object of the SPM class).
 */
class SPM
{
public:
    uint64_t nrows, ncols, nnz;
    
    SpmIterOrder type;
    
    SpmRowElem *elems__;
    uint64_t elems_size__;
    uint64_t *rowptr__;
    uint64_t rowptr_size__;
    
    bool elems_mapped;
    
    uint64_t GetNrRows()
    {
        return this->rowptr_size__ - 1;
    }

    SpmRowElem *RowBegin(uint64_t ridx=0);
    SpmRowElem *RowEnd(uint64_t ridx=0);

    uint64_t row_start;

    SPM() : type(NONE), elems__(NULL), rowptr__(NULL) {}
    ~SPM()
    {
        if (!elems_mapped && this->elems__)
            free(this->elems__);
        if (this->rowptr__)
            free(this->rowptr__);
    };

    class Builder;
    class VirtualBuilder;

    /**
     *  Function for filling the matrix using point iterators
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
     *  Load matrix from a file.
     *  @param filename name of the file that matrix is kept.
     *  @param in       buffer from which the matrix is taken.
     *  @return         spm class object with the characteristics of the matrix.
     */
    static SPM *LoadMMF(const char *filename);
    static SPM *LoadMMF(std::istream &in=std::cin);
    
    /**
     *  Load matrix from a file specifying the number of threads.
     *  @param filename name of the file that matrix is kept.
     *  @param in       buffer from which the matrix is taken.
     *  @param mmf      handler of MMF class.
     *  @param nr       number of threads to be used.
     *  @return         spm class object with the characteristics of the matrix.
     */
    static SPM *LoadMMF_mt(const char *mmf_file, const long nr);
    static SPM *LoadMMF_mt(std::istream &in, const long nr);
    static SPM *LoadMMF_mt(MMF &mmf, const long nr);

    void Print(std::ostream &out=std::cout);
    void PrintElems(std::ostream &out=std::cout);
    void PrintRows(std::ostream &out=std::cout);
    void PrintStats(std::ostream &out=std::cout);

    class PntIter;
    
    PntIter PointsBegin(uint64_t ridx=0);
    PntIter PointsEnd(uint64_t ridx=0);

    /**
     *  Transform matrix to specific type.
     *  @param type type of pattern to transform to.
     *  @param rs   line to start the transformation. If it is zero, it is the
     *              the beggining of matrix.
     *  @param re   line to end the transformation. If it is zero, it is the
     *              the end of matrix.
     */
    void Transform(SpmIterOrder type, uint64_t rs=0, uint64_t re=0);
    
    /**
     *  Find the reverse transform function. Transform from 'type' to 
     *  HORIZONTAL.
     *  @param type type of pattern to transform from.
     *  @return     reverse transform function.
     */
    TransformFn GetRevXformFn(SpmIterOrder type);
    
    /**
     *  Find the transform function. Transform from HORIZONTAL to 'type'.
     *  @param type type of pattern to transform to.
     *  @return     transform function.
     */
    TransformFn GetXformFn(SpmIterOrder type);
    
    /**
     *  Find the complete transform function. Tranform from 'from' to 'to'.
     *  @param from type of pattern to transform from.
     *  @param to   type of pattern to transform to.
     *  @return     complete transform function.
     */
    TransformFn GetTransformFn(SpmIterOrder from, SpmIterOrder to);
    
    /**
     *  Copy a window (sub-matrix) from matrix.
     *  @param rs     starting line of window.
     *  @param length number of lines to include in window.
     *  @return       spm object that includes the characteristics of the
     *                selected window.
     */
    SPM *ExtractWindow(uint64_t rs, uint64_t length);
    
    /**
     *  Get a window (sub-matrix) from matrix.
     *  @param rs     starting line of window.
     *  @param length number of lines to include in window.
     *  @return       spm object that includes the characteristics of the
     *                selected window.
     */
    SPM *GetWindow(uint64_t rs, uint64_t length);
    
    /**
     *  Copy the characteristics of the window (sub-matrix) back to to original
     *  matrix.
     *  @param window spm object that includes the characteristics of the
     *                window.
     */
    void PutWindow(const SPM *window);
};

void TestMMF(SPM *spm, const char *mmf_file);

/**
 *  Handles memory allocations during construction and destruction of spm 
 *  object.
 */
class SPM::Builder
{
public:
    SPM *spm;
    dynarray_t *da_elems;
    dynarray_t *da_rowptr;

    Builder(SPM *spm, uint64_t nr_elems=0, uint64_t nr_rows=0);
    virtual ~Builder();

    /**
     *  Allocates an element in matrix.
     *  @return pointer to element allocated.
     */
    virtual SpmRowElem *AllocElem();
    
    /**
     *  Allocates many elements in matrix.
     *  @param  nr number of elements allocated.
     *  @return    pointer to elements allocated.
     */
    virtual SpmRowElem *AllocElems(uint64_t nr);
    
    /**
     *  Counts elements allocated in matrix.
     *  @return number of elements allocated.
     */
    virtual uint64_t GetElemsCnt();
    
    /**
     *  Allocate elements for the next row.
     *  @param rdiff difference between the number of current and new line.
     */
    virtual void NewRow(uint64_t rdiff=1);
    
    virtual void Finalize();
};

/**
 *  Handles iterators that help in processing of the matrix.
 */
class SPM::PntIter : public std::iterator<std::forward_iterator_tag, CooElem>
{
public:
    SPM *spm;
    uint64_t row_idx;
    uint64_t elm_idx;
    
    PntIter();
    PntIter(SPM *s, uint64_t ridx);

    bool operator==(const PntIter &x);
    bool operator!=(const PntIter &x);
    void operator++();
    SpmCooElem operator*();
};

/**
 * Forward declarations for stream operations
 */
std::ostream &operator<<(std::ostream &os, const DeltaRLE::StatsVal &stats);
std::ostream &operator<<(std::ostream &os, const DeltaRLE &p);
std::ostream &operator<<(std::ostream &out, CooElem p);
std::ostream &operator<<(std::ostream &out, const SpmCooElem e);
std::ostream &operator<<(std::ostream &out, const SpmRowElem &elem);

}

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
