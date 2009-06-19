#ifndef CSX_SPM_H__
#define CSX_SPM_H__

#include <iostream>
#include <vector>
#include <iterator>

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

	class Generator;
	//virtual Generator generator(CooElem start) = 0;
};

class Pattern::Generator {
	public:
		virtual bool isEmpty() const = 0;
		virtual CooElem next() = 0;
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

} // csx namespace end

#endif /* CSX_SPM_H__ */
