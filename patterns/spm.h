#ifndef CSX_SPM_H__
#define CSX_SPM_H__

#include <vector>
#include <map>
#include <iterator>
#include <iostream>

#include <boost/function.hpp>

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

typedef enum {NONE=0, HORIZONTAL, VERTICAL, DIAGONAL, REV_DIAGONAL} SpmIterOrder;

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

	// Pattern detection and encoding
	class PrvData { };
	class Stats { };
	typedef std::map<PrvData *, Stats> StatsMap;
};

inline std::ostream &operator<<(std::ostream &os, const Pattern &p)
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
	void loadMMF(std::string mmf_file);
	void loadMMF(std::istream &in=std::cin);

	// Print Functions
	void Print(std::ostream &out=std::cout);
	void PrintRows(std::ostream &out=std::cout);

	template <typename IterT>
	void SetRows(IterT pnts_start, IterT pnts_end);

	// iterators that return a SpmCooElem
	class PointIter;
	PointIter points_begin();
	PointIter points_end();

	// same with PointIter, but removes elements
	class PointPoper;
	PointPoper points_pop_begin();
	PointPoper points_pop_end();

	// Transofrmation Functions
	TransformFn getRevXformFn(SpmIterOrder type);
	TransformFn getXformFn(SpmIterOrder type);
	TransformFn getTransformFn(SpmIterOrder from, SpmIterOrder to);
	void Transform(SpmIterOrder type);

	//
	static const long min_limit = 4;
	void DRLEncode();
	void DRLEncodeRow(SpmRowElems &oldrow, SpmRowElems &newrow);
	void doDRLEncode(uint64_t &col, std::vector<uint64_t> &xs, SpmRowElems &newrow);

	//
	void Draw(const char *filename, const int width=600, const int height=600);
};

} // csx namespace end

#endif /* CSX_SPM_H__ */
