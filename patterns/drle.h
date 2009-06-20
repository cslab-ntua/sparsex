#ifndef CSX_DRLE_H__
#define CSX_DRLE_H__

#include "spm.h"

#include <cassert>

namespace csx {

class DeltaRLE : public Pattern {
public:
	uint32_t size, delta;

	DeltaRLE(uint32_t _size, uint32_t _delta, SpmIterOrder _type):
	size(_size), delta(_delta){ ; }
	virtual DeltaRLE *clone() const
	{
		return new DeltaRLE(*this);
	}
	virtual long x_increase(SpmIterOrder order) const
	{
		long ret;
		ret = (order == this->type) ? (this->size*this->delta) : 1;
		return ret;
	}

	virtual std::ostream &print_on(std::ostream &out) const
	{
		out << "drle: size=" << this->size << " len=" << this->delta << " type=" << this->type;
		return out;
	}

	class Generator;
	Pattern::Generator *generator(CooElem start);
};

class DeltaRLE::Generator : public Pattern::Generator
{
	CooElem start;
	DeltaRLE *rle;
	long nr;

public:
	Generator(CooElem start_, DeltaRLE *rle_)
	: start(start_), rle(rle_), nr(0) { }

	virtual bool isEmpty() const {
		return (this->nr == this->rle->size);
	}
	virtual CooElem next() {
		CooElem ret(start);
		assert(this->nr <= this->rle->size);
		ret.x += (this->nr)*(this->rle->delta);
		this->nr += 1;
		return ret;
	}
};

Pattern::Generator *DeltaRLE::generator(CooElem start)
{
	DeltaRLE::Generator *g;
	g = new DeltaRLE::Generator(start, this);
	return g;
}

} // end of csx namespace

#endif /* CSX_DRLE_H__ */
