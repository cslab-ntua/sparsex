
#include "spm.h"
#include "drle.h"

#include <boost/foreach.hpp>
#define FOREACH BOOST_FOREACH

using namespace csx;

template <typename T>
struct RLE {
	long freq;
	T val;
};

template <typename T>
T DeltaEncode(T input)
{
	T output;
	typename T::iterator in, out;
	typename T::iterator::value_type prev, curr;

	output.resize(input.size());

	in = input.begin();
	out = output.begin();
	prev = *out++ = *in++;
	while (in < input.end()){
		curr = *in++;
		*out++ = curr - prev;
		prev = curr;
	}

	return output;
}

template <typename T>
std::vector< RLE<typename T::iterator::value_type> >
RLEncode(T input)
{
	typename T::iterator in;
	typename T::iterator::value_type curr;
	std::vector< RLE<typename T::iterator::value_type> > output;
	RLE<typename T::iterator::value_type> rle;

	in = input.begin();
	rle.freq = 1;
	rle.val = *in++;

	while (in < input.end()){
		curr = *in;
		if (rle.val == curr){
			rle.freq++;
		} else {
			output.push_back(rle);
			rle.freq = 1;
			rle.val = curr;
		}
		in++;
	}
	output.push_back(rle);
	return output;
}

void DRLE_Manager::updateStats(std::vector<uint64_t> &xs,
                               DeltaRLE::Stats &stats)
{
	std::vector< RLE<uint64_t> > rles;

	if (xs.size() == 0)
		return;

	rles = RLEncode(DeltaEncode(xs));
	FOREACH(RLE<uint64_t> &rle, rles){
		if (rle.freq >= this->min_limit){
			stats[rle.val].nnz += rle.freq;
			stats[rle.val].npatterns++;
		}
	}
	xs.clear();
}

DeltaRLE::Stats DRLE_Manager::generateStats()
{
	std::vector<uint64_t> xs;
	DeltaRLE::Stats stats;

	FOREACH(const SpmRowElems &row, this->spm.rows){
		FOREACH(const SpmRowElem &elem, row){
			if (elem.pattern == NULL){
				xs.push_back(elem.x);
				continue;
			}
			this->updateStats(xs, stats);
		}
		this->updateStats(xs, stats);
	}

	return stats;
}

void DRLE_Manager::doEncode(uint64_t &col,
                            std::vector<uint64_t> &xs,
                            SpmRowElems &newrow)
{
	std::vector< RLE<uint64_t> > rles;
	SpmRowElem elem;

	rles = RLEncode(DeltaEncode(xs));
	elem.pattern = NULL; // Default inserter (for push_back copies)
	FOREACH(RLE<uint64_t> &rle, rles){
		if (rle.freq >= this->min_limit){
			SpmRowElem *last_elem;

			col += rle.val; // go to the first
			elem.x = col;
			newrow.push_back(elem);
			last_elem = &newrow.back();
			last_elem->pattern = new DeltaRLE(rle.freq, rle.val, this->spm.type);
			last_elem = NULL;
			col += rle.val*(rle.freq - 1);
		} else {
			for (int i=0; i < rle.freq; i++){
				col += rle.val;
				elem.x = col;
				newrow.push_back(elem);
			}
		}
	}
}

void DRLE_Manager::EncodeRow(const SpmRowElems &oldrow,
                             SpmRowElems &newrow)
{
	std::vector<uint64_t> xs;
	uint64_t col;

	col = 0;
	FOREACH(const SpmRowElem &e, oldrow){
		if (e.pattern == NULL){
			xs.push_back(e.x);
			continue;
		}
		if (xs.size() != 0){
			doEncode(col, xs, newrow);
			xs.clear();
		}
		col += e.pattern->x_increase(this->spm.type);
		newrow.push_back(e);
	}
	if (xs.size() != 0){
		doEncode(col, xs, newrow);
		xs.clear();
	}
}

void DRLE_Manager::Encode()
{
	FOREACH(SpmRowElems &oldrow, this->spm.rows){
		SpmRowElems newrow;
		long newrow_size;
		// create new row
		EncodeRow(oldrow, newrow);
		// copy data
		newrow_size = newrow.size();
		oldrow.clear();
		oldrow.reserve(newrow_size);
		for (long i=0; i < newrow_size; i++){
			oldrow.push_back(newrow[i]);
		}
	}
}

Pattern::Generator *DeltaRLE::generator(CooElem start)
{
	DeltaRLE::Generator *g;
	g = new DeltaRLE::Generator(start, this);
	return g;
}

namespace csx {
void DRLE_OutStats(DeltaRLE::Stats &stats, SpmIdx &spm, std::ostream &os)
{
	DeltaRLE::Stats::iterator iter;
	for (iter=stats.begin(); iter != stats.end(); ++iter){
		os << "    " << iter->first << "-> "
		   << "np:" << iter->second.npatterns
		   << " nnz: " <<  100*((double)iter->second.nnz/(double)spm.nnz) << "%";
	}
}
} // end csx namespace

void DRLE_Manager::genAllStats()
{
	this->stats.clear();
	for (int t=HORIZONTAL; t != XFORM_MAX; t++){
		if (this->xforms[t])
			continue;

		SpmIterOrder type = SpmTypes[t];
		this->spm.Transform(type);
		this->stats[type] = this->generateStats();
		this->spm.Transform(HORIZONTAL);
	}
}
