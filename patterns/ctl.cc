#include <map>
#include <boost/foreach.h>

#define BOOST_FOREACH FOREACH

extern "C" {
	#include "dynarray.h"
}

#include "spm.h"
#include "ctl.h"

uint8_t CtlManager::getFlag(long pattern_id, uint64_t nnz)
{
	CtlFlags::PatMap::iterator pi;
	uint8_t ret;
	pi = this->patterns.find(pattern_id);
	if (pi == this->patterns.end()){
		ret = this->flag_avail++;
		CtlFlags::PatInfo patinfo(ret, nnz);
		this->patterns.insert(patinfo);
	} else {
		ret = pi->second.flag;
		pi->second.nr += nnz;
	}

	return ret;
}

// Static Pattern Id Mapping
//  8  -> delta 8
//  16 -> delta 16
//  32 -> delta 32
//  64 -> delta 64
//  10000 + delta => HORIZONTAL drle
//  20000 + delta => VERTICAL drle
//  30000 + delta => DIAGONAL drle
//  40000 + delta => REV_DIAGONAL drle
uint8_t *CtlManager::mkCtl()
{
	uint8_t *ret;
	dynarray_alloc(this->ctl_da, sizeof(uint8_t), 512);
	FOREACH(SpmRowElems &row, this->rows){
		this->doRow(row);
		row.clear();
	}

	ret = (uint8_t *)this->dynarray_destroy(this->ctl_da);
	return ret;
}

void CtlManager::doRow(const SpmRowElems &row)
{
	std::vector<uint64_t> xs;
	FOREACH(const SpmRowElem &spm_elem, row)
	{
		// check if this element contains a pattern
		if (spm_elem.pattern != NULL){
			if (xs.size() != 0)
				this->AddXs(xs)
			this->AddPattern(spm_elem.pattern);
			continue;
		}

		// check if we exceeded the maximum size for a unit
		assert(xs.size <= CTL_SIZE_MAX);
		if (xs.size() == CTL_SIZE_MAX){
			this->AddXs(xs);
			continue;
		}

		xs.push_back(row.x)
	}
}
