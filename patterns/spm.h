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

class MMF;

struct RowElem {
	uint64_t x;
	union {
		double val;
		double *vals;
	};
};

struct CooElem : public RowElem {
	uint64_t y;
};

static inline int CooCmp(const CooElem &p0, const CooElem &p1)
{
	int64_t ret;
	ret = p0.y - p1.y;
	if (ret == 0){
		ret = p0.x - p1.x;
	}
	if (ret > 0){
		return 1;
	} else if (ret < 0){
		return -1;
	} else {
		return 0;
	}
}

#define STRINGIFY__(s) #s
#define STRINGIFY(s)  STRINGIFY__(s)
#define BLOCK_ROW_TYPE_NAME(r)  BLOCK_R ## r
#define BLOCK_COL_TYPE_NAME(c)  BLOCK_C ## c


typedef enum {
	NONE=0,
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

//
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

/*
 *  Returns block alignment if `t' is a block type, otherwise returns `0'.
 */ 
static inline int isBlockType(SpmIterOrder t)
{
    /*
     *  This function assumes specific order for block type elements in
     *  SpmIterOrder enum (see the definition of SpmIterOrder).
     */
    if (t > BLOCK_TYPE_START && t < BLOCK_TYPE_END)
        return (t - BLOCK_COL_START > 0) ? (t - BLOCK_COL_START) :
            (t - BLOCK_TYPE_START);
    else
        return 0;
}

class Pattern {
public:
	SpmIterOrder type;

	virtual Pattern *clone() const = 0;
	virtual ~Pattern() {};
	virtual long getSize() const = 0;
	virtual long getPatId() const = 0;
	virtual long x_increase(SpmIterOrder spm_iter_order) const = 0;
	virtual uint64_t x_increase_jmp(SpmIterOrder spm_iter_order, uint64_t jmp) const = 0;
	virtual std::ostream &print_on(std::ostream &) const = 0;

	class Generator {
	public:
		virtual bool isEmpty() const = 0;
		virtual CooElem next() = 0;
	};
	virtual Generator *generator(CooElem start) = 0;

	// get the next x
	// Todo: maybe use this for a common generator
	virtual uint64_t getNextX(uint64_t x0) const = 0;

	// stats for a specific pattern. For now it's just the number of non-zero
	// elements that adhere to this pattern.
	class StatsVal {
	public:
		uint64_t nnz;
		long npatterns;
		//StatsVal() : nnz(0) { }
	};
};


class SpmPattern {
public:
	Pattern *pattern;

	SpmPattern(void) : pattern(NULL) { ; }
	SpmPattern(const SpmPattern &spm_p){
		// make a copy, so we don't fsck up when destructor is called
		Pattern *p;
		this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->clone();
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
		this->pattern = ((p = spm_p.pattern) == NULL) ? NULL : p->clone();
		return *this;
	}
};

class SpmCooElem: public CooElem, public SpmPattern {};
class SpmRowElem: public RowElem, public SpmPattern {};

void mk_row_elem(const CooElem &p, SpmRowElem *ret);
void mk_row_elem(const SpmCooElem &p, SpmRowElem *ret);
void mk_row_elem(const SpmRowElem &p, SpmRowElem *ret);

#if 0
typedef std::vector<CooElem> SpmPoints;
typedef std::vector<SpmCooElem> SpmCooElems;
typedef std::vector<SpmRowElem> SpmRowElems;
typedef std::vector<SpmRowElems> SpmRows;
#endif

typedef std::iterator<std::forward_iterator_tag, CooElem> SpmPointIter;
typedef boost::function<void (CooElem &p)> TransformFn;


class SPM {
public:
	// Charateristics of sparse matrix
	uint64_t nrows, ncols, nnz;

	SpmIterOrder type;

	// Elements of sparse matrix
	SpmRowElem *elems__;
	uint64_t elems_size__;
	uint64_t *rowptr__;
	uint64_t rowptr_size__;

	uint64_t getNrRows() { return this->rowptr_size__ - 1; }
	//// SpmRowElem iterators
	SpmRowElem *rbegin(uint64_t ridx=0); // start of ridx row
	SpmRowElem *rend(uint64_t ridx=0); // end of ridx row

	// Since this can be a partition of the original matrix,
	// this designates the first row of this partition
	uint64_t row_start;

	SPM() : type(NONE), elems__(NULL), rowptr__(NULL) {}
	~SPM(){
		if (this->elems__)
			free(this->elems__);

		if (this->rowptr__)
			free(this->rowptr__);
	};

	// low level creation of spm
	class Builder;

	// function for filling the matrix using point iterators
	// pnts_start, pnts_end : point iterators start/end
	// first_row: first row of the matrix (starting from 1)
	// limit : limit for the number of elements. When this limit is reached the
	//         function returns with the number of elements used.
	// elems_nr, rows_nr : size for elems and rows (initial) allocation
	template <typename IterT>
	uint64_t SetElems(IterT &pnts_start, const IterT &pnts_end, uint64_t first_row,
	                  unsigned long limit=0,
	                  uint64_t elems_nr=0,
	                  uint64_t rows_nr=0);

	// load matrix from an MMF file
	static SPM *loadMMF(const char *filename);
	static SPM *loadMMF(std::istream &in=std::cin);
	static SPM *loadMMF_mt(const char *mmf_file, const long nr);
	static SPM *loadMMF_mt(std::istream &in, const long nr);
	static SPM *loadMMF_mt(MMF &mmf, const long nr);

	// Print Functions
	void Print(std::ostream &out=std::cout);
	void PrintElems(std::ostream &out=std::cout);
	void PrintRows(std::ostream &out=std::cout);
    	void PrintStats(std::ostream &out=std::cout);

	// iterators of the sparse matrix that return a SpmCooElem
	class PntIter;
	PntIter points_begin(uint64_t ridx=0);
	PntIter points_end(uint64_t ridx=0);

	// Transformation Functions
	void Transform(SpmIterOrder type, uint64_t rs=0, uint64_t re=0);
	TransformFn getRevXformFn(SpmIterOrder type);
	TransformFn getXformFn(SpmIterOrder type);
	TransformFn getTransformFn(SpmIterOrder from, SpmIterOrder to);
};

void TestMMF(SPM *spm, const char *mmf_file);

class SPM::Builder
{
public:
	SPM *spm;
	dynarray_t *da_elems;
	dynarray_t *da_rowptr;

	Builder(SPM *spm, uint64_t elems_nr=0, uint64_t rows_nr=0);
	~Builder();

	// Using these two functions requires caution: Memory is uninitialized
	// use mk_row_elem functions
	SpmRowElem *AllocElem();
	SpmRowElem *AllocElems(uint64_t nr);

	uint64_t getElemsCnt();
	void newRow(uint64_t rdiff=1);
	void Finalize();
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
#if 0
std::ostream &operator<<(std::ostream &out, const SpmRowElems &elems);
std::ostream &operator<<(std::ostream &out, const SpmRows &rows);
std::ostream &operator<<(std::ostream &out, SpmIdx::PointIter pi);
#endif

} // csx namespace end

#endif /* CSX_SPM_H__ */
