#ifndef CSX_SPM_H__
#define CSX_SPM_H__

#include <vector>
#include <map>
#include <iterator>
#include <iostream>
#include <iomanip>

#include <boost/function.hpp>
#include <boost/foreach.hpp>

namespace csx {

class CooElem {
public:
	uint64_t y;
	uint64_t x;
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

class RowElem {
public:
	uint64_t x;
};

typedef enum {
	NONE=0,
	HORIZONTAL,
	VERTICAL,
	DIAGONAL,
	REV_DIAGONAL,
	XFORM_MAX
} SpmIterOrder;

// !#@$$#!@#%!$#%$!$%!!$%$%
const SpmIterOrder SpmTypes[] = {
	NONE,
	HORIZONTAL,
	VERTICAL,
	DIAGONAL,
	REV_DIAGONAL,
	XFORM_MAX
};

const char *SpmTypesNames[] = {
	"__NONE__",
	"HORIZONTAL",
	"VERTICAL",
	"DIAGONAL",
	"REV_DIAGONAL",
	"__XFORM_MAX__"
};

class Pattern {
public:
	SpmIterOrder type;

	virtual Pattern *clone() const = 0;
	//virtual ~Pattern() const = 0;
	virtual long x_increase(SpmIterOrder spm_iter_order) const = 0;
	virtual std::ostream &print_on(std::ostream &) const = 0;

	class Generator {
	public:
		virtual bool isEmpty() const = 0;
		virtual CooElem next() = 0;
	};
	virtual Generator *generator(CooElem start) = 0;

	// stats for a specific pattern. For now it's just the number of non-zero
	// elements that adhere to this pattern.
	class StatsVal {
	public:
		uint64_t nnz;
		long npatterns;
		//StatsVal() : nnz(0) { }
	};
};

static inline std::ostream &operator<<(std::ostream &os, const Pattern::StatsVal &stats)
{
	os << "nnz: " << stats.nnz;
	return os;
}

static inline std::ostream &operator<<(std::ostream &os, const Pattern &p)
{
	os << " (";
	p.print_on(os);
	os << " type:" << p.type << ") ";
	return os;
}

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

typedef std::vector<CooElem> SpmPoints;
typedef std::iterator<std::forward_iterator_tag, CooElem> SpmPointIter;
typedef std::vector<SpmCooElem> SpmCooElems;
typedef std::vector<SpmRowElem> SpmRowElems;
typedef std::vector<SpmRowElems> SpmRows;
typedef boost::function<void (CooElem &p)> TransformFn;

class SpmIdx {
public:
	uint64_t nrows, ncols, nnz;
	SpmIterOrder type;
	SpmRows rows;

	SpmIdx() {type = NONE;};
	~SpmIdx() {};

	// load matrix from an MMF file
	void loadMMF(const char *filename);
	void loadMMF(std::istream &in=std::cin);

	// Print Functions
	void Print(std::ostream &out=std::cout);
	void PrintRows(std::ostream &out=std::cout);

	template <typename IterT>
	void SetRows(IterT pnts_start, IterT pnts_end);

	// iterators that return a SpmCooElem
	class PointIter;
	PointIter points_begin(uint64_t ridx=0, uint64_t eidx=0);
	PointIter points_end(uint64_t ridx=0, uint64_t eidx=0);

	// same with PointIter, but removes elements
	class PointPoper;
	PointPoper points_pop_begin(uint64_t ridx=0, uint64_t eidx=0);
	PointPoper points_pop_end(uint64_t ridx=0, uint64_t eidx=0);

	typedef SpmRows::iterator RowIter;
	RowIter rows_iter(uint64_t r){
		RowIter ret;
		ret = this->rows.begin();
		std::advance(ret, r);
		return ret;
	}
	RowIter rbegin() { return  rows.begin(); }
	RowIter rend() { return rows.end(); }

	// Transofrmation Functions
	TransformFn getRevXformFn(SpmIterOrder type);
	TransformFn getXformFn(SpmIterOrder type);
	TransformFn getTransformFn(SpmIterOrder from, SpmIterOrder to);
	void Transform(SpmIterOrder type, uint64_t rs=0, uint64_t re=0);
};

class SpmIdx::PointIter
: public std::iterator<std::forward_iterator_tag, CooElem>
{
public:
	uint64_t row_idx;
	uint64_t elm_idx;
	SpmIdx   *spm;
	PointIter(): row_idx(0), elm_idx(0), spm(NULL) {;}
	PointIter(uint64_t ridx, uint64_t eidx, SpmIdx *s)
	: row_idx(ridx), elm_idx(eidx), spm(s) {
		// find the first non zero row
		while (row_idx < spm->nrows && spm->rows[row_idx].size() == 0) {
			row_idx++;
		}
	}

	bool operator==(const PointIter &x){
		return (spm == x.spm)
			   && (row_idx == x.row_idx)
			   && (elm_idx == x.elm_idx);
	}

	bool operator!=(const PointIter &x){
		return !(*this == x);
	}

	void operator++() {
		uint64_t rows_nr, row_elems_nr;

		rows_nr = spm->rows.size();
		row_elems_nr = spm->rows[row_idx].size();

		// equality means somebody did a ++ on an ended iterator
		assert(row_idx < rows_nr);

		elm_idx++;
		assert(elm_idx <= row_elems_nr);
		if (elm_idx < row_elems_nr){
			return ;
		}

		// change row
		do {
			elm_idx = 0;
			row_idx++;
		} while (row_idx < spm->rows.size() && spm->rows[row_idx].size() == 0);
	}

	SpmCooElem operator*(){
		SpmCooElem ret;
		ret.y = row_idx + 1;
		ret.x = spm->rows[row_idx][elm_idx].x;
		Pattern *p = spm->rows[row_idx][elm_idx].pattern;
		ret.pattern = (p == NULL) ? NULL : p->clone();
		return ret;
	}
};


class SpmIdx::PointPoper : public SpmIdx::PointIter
{
public:
	PointPoper() : PointIter() { ; }
	PointPoper(uint64_t ridx, uint64_t eidx, SpmIdx *s)
	: PointIter(ridx, eidx, s) { ; }

	void operator++() {
		uint64_t rows_nr, row_elems_nr;

		rows_nr = spm->rows.size();
		row_elems_nr = spm->rows[row_idx].size();

		// equality means somebody did a ++ on an ended iterator
		assert(row_idx < rows_nr);

		elm_idx++;
		assert(elm_idx <= row_elems_nr);
		if (elm_idx < row_elems_nr){
			return;
		}

		// remove all elements from previous row
		spm->rows[row_idx].resize(0);
		// change row
		do {
			elm_idx = 0;
			row_idx++;
		} while (row_idx < rows_nr && spm->rows[row_idx].size() == 0);
	}
};


} // csx namespace end

#endif /* CSX_SPM_H__ */
