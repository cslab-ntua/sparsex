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

namespace csx {

class MMF;  //*** Create a file with possible forward declarations?

#define STRINGIFY__(s) #s
#define STRINGIFY(s)  STRINGIFY__(s)
#define BLOCK_ROW_TYPE_NAME(r)  BLOCK_R ## r
#define BLOCK_COL_TYPE_NAME(c)  BLOCK_C ## c

/**
 *  SpmIterOrder enum.
 *  Order patterns used in CSX and match them with numbers.
 */
typedef enum {
    NONE = 0,
    HORIZONTAL,
    VERTICAL,
    DIAGONAL,
    REV_DIAGONAL,
    BLOCK_TYPE_START,           //*** possible insertion of BLOCK_ROW_START
    BLOCK_ROW_TYPE_NAME(1),     //*** in order to be symmetric.
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

const SpmIterOrder SpmTypes[] = {
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

const char *SpmTypesNames[] = {
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
 * IsBlockType function.
 * @param t a type of pattern.
 * @see     enum SpmIterOrder.
 * @return  block alignment if 't' is a block type,
 *          0 otherwise.
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

//*** Maybe is better to do RowElem and CooElem classes?
/**
 *  RowElem struct.
 *  Holds column and value or values of CSX elements.
 */
struct RowElem {
    uint64_t x;
    union {
        double val;
        double *vals;
    };
};

/**
 *  CooElem struct.
 *  Expansion of RowElem which holds row of CSX elements.
 */
struct CooElem : public RowElem {
    uint64_t y;
};

/**
 * CooCmp function.
 * @param p0 an element with row, column and value or values.
 * @param p1 an element with row, column and value or values.
 * @see      struct CooElem.
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
 *  Pattern Class.
 *  Class that describes pattern of non-zero elements.
 */
class Pattern {
public:
    SpmIterOrder type;

    /**
     *  A virtual Clone function.
     *  @return a copy of the chosen object of class Pattern.
     */
    virtual Pattern *Clone() const = 0;
    
    /**
     *  A destructor of class Pattern.
     */
    virtual ~Pattern() {};
    
    /**
     *  A virtual GetSize function.
     *  @return the size of the selected object of class Pattern.
     */
    virtual long GetSize() const = 0;
    
    /**
     *  A virtual GetPatId function.
     *  @return a unique id of the type of the selected object of class Pattern.
     */
    virtual long GetPatId() const = 0;
    
    //*** We do not use ColIncrease.
    /**
     *  A virtual ColIncrease function.
     *  @param spm_iter_order a type of pattern.
     *  @return               number of columns to procceed in order to go to
     *                        the next CSX element.
     */
    virtual long ColIncrease(SpmIterOrder spm_iter_order) const = 0;
    
    /**
     *  A virtual ColIncrease function.
     *  @param spm_iter_order a type of pattern.
     *  @param jmp            a number of columns to jump.
     *  @return               number of columns to procceed in order to go to
     *                        the next CSX element.
     */
    virtual uint64_t ColIncreaseJmp(SpmIterOrder spm_iter_order, uint64_t jmp) 
                                    const = 0;

    virtual std::ostream &PrintOn(std::ostream &) const = 0;

    //*** We use class Generator only on PrintElems which is not used.
    class Generator {
    public:
        virtual bool isEmpty() const = 0;
        virtual CooElem next() = 0;
    };
    virtual Generator *generator(CooElem start) = 0;

    /**
     *  Get the column of the next non-zero element.
     */
    //*** I do not understand todo operation.
    // TODO: maybe use this for a common generator
    virtual uint64_t GetNextCol(uint64_t x0) const = 0;
    
    /**
     *  Stats for a specific type of pattern of sparse matrix. More
     *  specifically, it just holds the number of non-zero elements and 
     *  patterns that adhere to this specific type of pattern.
     */
    class StatsVal {
    public:
        uint64_t nnz;
        long npatterns;
        
        void Update(const StatsVal& new_vals) {
            this->nnz += new_vals.nnz;
            this->npatterns += new_vals.npatterns;
        }
    };
};

//*** Merge SpmPattern with Pattern Class? 
class SpmPattern
{
public:
    Pattern *pattern;
    
    SpmPattern(void) : pattern(NULL) { }
    SpmPattern(const SpmPattern &spm_p){
        // make a copy, so we don't fsck up when destructor is called
        Pattern *p;
        this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->Clone();
    }
    
    ~SpmPattern() {
        if (this->pattern != NULL)
            delete this->pattern;
    }
    
    SpmPattern& operator=(const SpmPattern &spm_p) {
        if (this == &spm_p)
            return *this;
            
        if (this->pattern != NULL)
            delete this->pattern;
            
        // make a copy, so we don't fsck up when destructor is called
        Pattern *p;
        this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->Clone();
        return *this;
    }
};

/**
 *  Expansion of CooElem and RowElem which holds pattern info about the CSX
 *  elements.
 */
class SpmCooElem: public CooElem, public SpmPattern {};
class SpmRowElem: public RowElem, public SpmPattern {};

/**
 *  Copies element 'p' to SpmRowElem 'ret'.
 */
//*** This functions are not used (the first and the second one).
void MakeRowElem(const CooElem &p, SpmRowElem *ret);
void MakeRowElem(const SpmCooElem &p, SpmRowElem *ret);

void MakeRowElem(const SpmRowElem &p, SpmRowElem *ret);

//*** The SpmPointIter type is not used.
typedef std::iterator<std::forward_iterator_tag, CooElem> SpmPointIter;
typedef boost::function<void (CooElem &p)> TransformFn;

class SPM
{
public:
    // Charateristics of sparse matrix
    uint64_t nrows, ncols, nnz;
    SpmIterOrder type;
    
    // Elements of sparse matrix
    SpmRowElem *elems__;
    uint64_t elems_size__;
    uint64_t *rowptr__;
    uint64_t rowptr_size__;
    bool elems_mapped;
    
    uint64_t GetNrRows()
    {
        return this->rowptr_size__ - 1;
    }
    
    // SpmRowElem iterators
    SpmRowElem *RowBegin(uint64_t ridx=0); 	// start of ridx row
    SpmRowElem *RowEnd(uint64_t ridx=0); 	// end of ridx row
    
    // Since this can be a partition of the original matrix,
    // this designates the first row of this partition
    uint64_t row_start;
    
    SPM() : type(NONE), elems__(NULL), rowptr__(NULL) {}
    ~SPM()
    {
        if (!elems_mapped && this->elems__)
            free(this->elems__);
        if (this->rowptr__)
            free(this->rowptr__);
    };
    
    // low level creation of spm
    class Builder;
    class VirtualBuilder;
    
    // function for filling the matrix using point iterators
    // pnts_start, pnts_end : point iterators start/end
    // first_row: first row of the matrix (starting from 1)
    // limit : limit for the number of elements. When this limit is reached the
    //         function returns with the number of elements used.
    // elems_nr, rows_nr : size for elems and rows (initial) allocation
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit = 0,
                      uint64_t elems_nr = 0, uint64_t rows_nr = 0);
                      
    template <typename IterT>
    uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end,
                      uint64_t first_row, unsigned long limit,
                      uint64_t elems_nr, uint64_t rows_nr,
                      SPM::Builder *SpmBld);
                      
    // Load matrix from an MMF file
    static SPM *LoadMMF(const char *filename);
    static SPM *LoadMMF(std::istream &in=std::cin);
    static SPM *LoadMMF_mt(const char *mmf_file, const long nr);
    static SPM *LoadMMF_mt(std::istream &in, const long nr);
    static SPM *LoadMMF_mt(MMF &mmf, const long nr);
    
    // Print Functions
    //*** Print functions are not used.
    void Print(std::ostream &out=std::cout);
    void PrintElems(std::ostream &out=std::cout);
    void PrintRows(std::ostream &out=std::cout);
    void PrintStats(std::ostream &out=std::cout);
    
    // iterators of the sparse matrix that return a SpmCooElem
    class PntIter;
    PntIter PointsBegin(uint64_t ridx=0);
    PntIter PointsEnd(uint64_t ridx=0);
    
    // Transformation Functions
    void Transform(SpmIterOrder type, uint64_t rs=0, uint64_t re=0);
    TransformFn GetRevXformFn(SpmIterOrder type);
    TransformFn GetXformFn(SpmIterOrder type);
    TransformFn GetTransformFn(SpmIterOrder from, SpmIterOrder to);
    
    //*** ExtractWindow function is not used.
    SPM *ExtractWindow(uint64_t rs, uint64_t length);
    SPM *GetWindow(uint64_t rs, uint64_t length);
    void PutWindow(const SPM *window);
};

void TestMMF(SPM *spm, const char *mmf_file);

class SPM::Builder
{
public:
    SPM *spm;
    dynarray_t *da_elems;
    dynarray_t *da_rowptr;
    
    Builder(SPM *spm, uint64_t elems_nr=0, uint64_t rows_nr=0);
    virtual ~Builder();
    
    // Using these two functions requires caution: Memory is uninitialized
    // use MakeRowElem() functions
    virtual SpmRowElem *AllocElem();
    virtual SpmRowElem *AllocElems(uint64_t nr);
    
    virtual uint64_t GetElemsCnt();
    virtual void NewRow(uint64_t rdiff=1);
    virtual void Finalize();
};

class SPM::PntIter : public std::iterator<std::forward_iterator_tag, CooElem>
{
public:
    SPM *spm;
    uint64_t row_idx; // current row
    uint64_t elm_idx; // current element
    
    // Empty initialization
    PntIter();
    // Initialization for a matrix and start at the specified row
    PntIter(SPM *s, uint64_t ridx);
    
    bool operator==(const PntIter &x);
    bool operator!=(const PntIter &x);
    void operator++();
    SpmCooElem operator*();
};

// Forward declarations for stream operations
std::ostream &operator<<(std::ostream &os, const Pattern::StatsVal &stats);
std::ostream &operator<<(std::ostream &os, const Pattern &p);
std::ostream &operator<<(std::ostream &out, CooElem p);
std::ostream &operator<<(std::ostream &out, const SpmCooElem e);
std::ostream &operator<<(std::ostream &out, const SpmRowElem &elem);

} // csx namespace end

#endif /* CSX_SPM_H__ */
