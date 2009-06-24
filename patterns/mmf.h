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
	bool valid;

public:
	void _set(){
		double v;
		this->valid = this->mmf->next(this->elem.y, this->elem.x, v);
	}
	iterator() { assert(false); }
	iterator(MMF *mmf_, uint64_t cnt_): mmf(mmf_), cnt(cnt_)
	{
		// this is the initializer
		if (this->cnt == 0){
			this->_set();
		}
	}


	bool operator==(const iterator &i){
		//std::cout << "me: " << this->mmf << " " << this->cnt
		//          << " i: " << i.mmf << " " << i.cnt << "\n";
		return (this->mmf == i.mmf) && (this->cnt == i.cnt);
	}

	bool operator!=(const iterator &i){
		return !(*this == i);
	}

	void operator++(){
		this->cnt++;
		this->_set();
	}

	SpmCooElem operator*(){
		if (!this->valid){
			std::cout << "Requesting dereference, but mmf ended\n"
			          << "cnt: " << this->cnt << std::endl;
			assert(false);
		}
		assert(valid);
		return this->elem;
	}
};

MMF::iterator MMF::begin() { return MMF::iterator(this, 0); }
MMF::iterator MMF::end() { return MMF::iterator(this, this->nnz); }

} // csx namespace end

#endif
