#ifndef MMF_H__
#define MMF_H__

#include "spm.h"

#include <iostream>
#include <iterator>

namespace csx {

class MMF
{
private:
	std::istream &in;

public:
	uint64_t nrows, ncols, nnz;

	// initialization
	MMF(std::istream &in);

	// get next element (false if end)
	bool next(uint64_t &y, uint64_t &x, double &val);

	// CooElem iterator
	class iterator;
	iterator begin();
	iterator end();
};

class MMF::iterator :
public std::iterator<std::forward_iterator_tag, CooElem>
{
	MMF *mmf;
	uint64_t cnt;
	SpmCooElem elem;

public:
	iterator(MMF *mmf_, uint64_t cnt_): mmf(mmf_), cnt(cnt_) {}

	bool operator==(const iterator &i)
	{
		return (this->mmf == i.mmf) && (this->cnt == i.cnt);
	}

	bool operator!=(const iterator &i)
	{
		return !(*this == i);
	}

	void operator++()
	{
		bool ret;
		double v;
		ret = this->mmf->next(this->elem.y, this->elem.x, v);
		assert(ret);
		this->cnt++;
	}

	SpmCooElem operator*()
	{
		return this->elem;
	}
};

MMF::iterator MMF::begin() { return MMF::iterator(this, 0); }
MMF::iterator MMF::end() { return MMF::iterator(this, this->nnz); }

} // csx namespace end

#endif
