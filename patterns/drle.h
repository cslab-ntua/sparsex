#ifndef CSX_DRLE_H__
#define CSX_DRLE_H__

#include <map>
#include <bitset>

#include "spm.h"

#include <cassert>

namespace csx {

class DeltaRLE : public Pattern {
	uint32_t size, delta;

public:
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

	virtual uint64_t x_increase_jmp(SpmIterOrder order, uint64_t jmp) const
	{
		long ret;

		ret = jmp;
		if (order == this->type){
			ret += (this->size*this->delta);
		}

		return ret;
	}

	virtual std::ostream &print_on(std::ostream &out) const
	{
		out << "drle: size=" << this->size << " len=" << this->delta << " type=" << this->type;
		return out;
	}

	virtual long getPatId() const
	{
		//******** Static Pattern Id Mapping:
		//  10000 + delta => HORIZONTAL drle
		#define PID_HORIZ_BASE 10000
		//  20000 + delta => VERTICAL drle
		#define PID_VERT_BASE 20000
		//  30000 + delta => DIAGONAL drle
		#define PID_DIAG_BASE 30000
		//  40000 + delta => REV_DIAGONAL drle
		#define PID_rDIAG_BASE 40000
		switch (type){
			case HORIZONTAL: return PID_HORIZ_BASE + this->delta;
			case VERTICAL: return PID_VERT_BASE + this->delta;
			case DIAGONAL: return PID_DIAG_BASE + this->delta;
			case REV_DIAGONAL: return PID_rDIAG_BASE + this->delta;
			default: assert(false);
		}
	}

	virtual long getSize() const
	{
		return this->size;
	}

	class Generator;
	Pattern::Generator *generator(CooElem start);

	// key => delta value of rle
	typedef std::map<uint64_t, Pattern::StatsVal> Stats;
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

class DRLE_Manager {
public:
	SpmIdx *spm;
	long min_limit; // minimum length for RLEs
	long max_limit; // maximum length for RLEs
	double min_perc; // min nnz percentage for considering an RLE

	DRLE_Manager(SpmIdx *_spm,
	             long min_limit_=4, long max_limit_ = 255,
	             double min_perc_=.1)
	: spm(_spm), min_limit(min_limit_), max_limit(max_limit_), min_perc(min_perc_) { }

	DeltaRLE::Stats generateStats();

	typedef std::map <SpmIterOrder, DeltaRLE::Stats> StatsMap;
	StatsMap stats;
	void genAllStats();
	void outStats(std::ostream &os=std::cout);

	std::bitset<XFORM_MAX> xforms;
	void Encode();

	SpmIterOrder chooseType();
	uint64_t getTypeNNZ(SpmIterOrder type);

	// create the ctl array
	// (this will destroy the spm)
	unsigned char *mkCtl();

private:
	void doEncode(uint64_t &col, std::vector<uint64_t> &xs, SpmRowElems &newrow);
	void EncodeRow(const SpmRowElems &oldrow, SpmRowElems &newrow);
	void updateStats(std::vector<uint64_t> &xs, DeltaRLE::Stats &stats);
};

#if 0
inline std::ostream &operator<<(std::ostream &os, const DeltaRLE::Stats &stats)
{
	DeltaRLE::Stats::const_iterator iter;
	for (iter=stats.begin(); iter != stats.end(); ++iter){
		os << iter->first << " " << iter->second << "  ";
	}
	return os;
}
#endif

void DRLE_OutStats(DeltaRLE::Stats &stats, SpmIdx &spm, std::ostream &os);

} // end of csx namespace

#endif /* CSX_DRLE_H__ */
