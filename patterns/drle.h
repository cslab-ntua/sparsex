#ifndef CSX_DRLE_H__
#define CSX_DRLE_H__

#include <map>
#include <set>
#include <bitset>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <stdexcept>

#include "spm.h"

namespace csx {

class DeltaRLE : public Pattern {
	uint32_t size, delta;

public:
	DeltaRLE(uint32_t _size, uint32_t _delta, SpmIterOrder _type):
	size(_size), delta(_delta){ this->type = _type; }

	virtual DeltaRLE *clone() const
	{
		return new DeltaRLE(*this);
	}

	virtual long x_increase(SpmIterOrder order) const
	{
		long        ret;

		ret = (order == this->type) ? (this->size*this->delta) : 1;
		return ret;
	}

	virtual uint64_t x_increase_jmp(SpmIterOrder order, uint64_t jmp) const
	{
		long ret;

		ret = jmp;
		if (order == this->type){
			ret += ((this->size-1)*this->delta);
		}

		return ret;
	}

	virtual std::ostream &print_on(std::ostream &out) const
	{
		out << "drle: size=" << this->size << " len=" << this->delta << " type=" << this->type;
		return out;
	}

	#define PID_OFFSET 10000
	virtual long getPatId() const
	{
		assert(this->type > NONE && this->type < BLOCK_TYPE_START);
		return type*PID_OFFSET + this->delta;

// 		//******** Static Pattern Id Mapping:
// 		//  10000 + delta => HORIZONTAL drle
// 		#define PID_HORIZ_BASE 10000
// 		//  20000 + delta => VERTICAL drle
// 		#define PID_VERT_BASE 20000
// 		//  30000 + delta => DIAGONAL drle
// 		#define PID_DIAG_BASE 30000
// 		//  40000 + delta => REV_DIAGONAL drle
// 		#define PID_rDIAG_BASE 40000
// 		switch (type){
// 			case HORIZONTAL: return PID_HORIZ_BASE + this->delta;
// 			case VERTICAL: return PID_VERT_BASE + this->delta;
// 			case DIAGONAL: return PID_DIAG_BASE + this->delta;
// 			case REV_DIAGONAL: return PID_rDIAG_BASE + this->delta;
// 			default: assert(false);
// 		}
	}

	virtual long getSize() const
	{
		return this->size;
	}

	virtual uint64_t getNextX(uint64_t x0) const
	{
        	return (x0 + this->delta);
	}

	class Generator;
	Pattern::Generator *generator(CooElem start);
	// key => delta value of rle
	// key => second dimension for block rles.
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
		ret.x += (this->nr)*this->rle->delta;
		this->nr += 1;
		return ret;
	}

};

class DRLE_Manager {
public:
	SPM *spm;
	long min_limit; // minimum length for RLEs
	long max_limit; // maximum length for RLEs
	double min_perc; // min nnz percentage for considering an RLE

	typedef enum {
		SPLIT_BY_ROWS = 0,
		SPLIT_BY_NNZ,
	} split_alg_t;

        DRLE_Manager(SPM *_spm,
             long min_limit_=4,
             long max_limit_ = std::numeric_limits<long>::max(),
             double min_perc_=.1, uint64_t sort_window_size_ = 0,
             split_alg_t split_type_ = SPLIT_BY_ROWS,
             double probability = 1.0);

	DeltaRLE::Stats generateStats(SPM *spm, uint64_t rs, uint64_t re);
	DeltaRLE::Stats generateStats(uint64_t rs, uint64_t re);

	typedef std::map <SpmIterOrder, DeltaRLE::Stats> StatsMap;
	StatsMap stats;
	void genAllStats();
	void outStats(std::ostream &os=std::cout);
	void outStats(char *buffer);
    
	std::map <SpmIterOrder, std::set<uint64_t> > DeltasToEncode;
	std::bitset<XFORM_MAX> xforms_ignore;

	void addIgnore(SpmIterOrder type);
	void ignoreAll();
	void removeIgnore(SpmIterOrder type);
	void removeAll();

	SpmIterOrder chooseType();
	uint64_t getTypeScore(SpmIterOrder type);

	void Encode(SpmIterOrder type=NONE);
	void Decode(SpmIterOrder type=NONE);			//Prosthesa edw
	void EncodeAll(char *buffer);
	void MakeEncodeTree();					//edw
	void EncodeSerial(int *xform_str);			//edw
    void set_sampling_probability(double probability) {
        if (probability < 0.0 || probability > 1.0)
            throw new std::invalid_argument("invalid sampling probability");

        sampling_probability = probability;
    }

    double get_sampling_probability(double probability) {
        return sampling_probability;
    }

private:
	void doEncode(std::vector<uint64_t> &xs,
	              std::vector<double> &vs,
	              std::vector<SpmRowElem> &newrow);

	void doEncodeBlock(std::vector<uint64_t> &xs,
                       std::vector<double> &vs,
                       std::vector<SpmRowElem> &newrow);
				
	void doDecode(const SpmRowElem *elem,			//edw
		      std::vector<SpmRowElem> &newrow);

	void EncodeRow(const SpmRowElem *rstart,
	               const SpmRowElem *rend,
	               std::vector<SpmRowElem> &newrow);


	void DecodeRow(const SpmRowElem *rstart,
		       const SpmRowElem *rend, 
		       std::vector<SpmRowElem> &newrow);	//kai edw


	void updateStats(SPM *spm, std::vector<uint64_t> &xs,
			 DeltaRLE::Stats &stats);
	
	void updateStats2(SPM *spm, std::vector<uint64_t> &xs,
			 DeltaRLE::Stats &stats);

	void updateStats(std::vector<uint64_t> &xs,
	                 DeltaRLE::Stats &stats);
    	
	void updateStatsBlock(std::vector<uint64_t> &xs,
                          DeltaRLE::Stats &stats, uint64_t align);

    void updateStats(SpmIterOrder type, DeltaRLE::Stats stats);
    void correct_stats(SpmIterOrder type, double factor);

    void compute_sort_splits();
    void check_and_set_sorting();
    void do_compute_sort_splits_by_rows();
    void do_compute_sort_splits_by_nnz();
    void do_check_sort_by_rows();
    void do_check_sort_by_nnz();

    bool sort_windows;
    uint64_t sort_window_size;
    split_alg_t split_type;
    std::vector<uint64_t> sort_splits;
    typedef std::vector<uint64_t>::iterator sort_split_iterator;
    double sampling_probability;
};

class BlockRLE : public DeltaRLE {
	uint32_t    other_dim;

public:
	BlockRLE(uint32_t size_, uint32_t other_dim_, SpmIterOrder type_)
	: DeltaRLE(size_, 1, (assert(isBlockType(type_)), type_))
	{
		this->other_dim = other_dim_;
        }

	virtual uint32_t getOtherDim() const
	{
		return this->other_dim;
	}

	virtual long getPatId() const
	{
		return PID_OFFSET*this->type + this->other_dim;
	}

	virtual BlockRLE *clone() const
	{
		return new BlockRLE(*this);
	}
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

void DRLE_OutStats(DeltaRLE::Stats &stats, SPM &spm, std::ostream &os);
void DRLE_OutStats(DeltaRLE::Stats &stats, SPM &spm, char *buffer);

class Node {								//Prosthesa edw
public:
	uint32_t depth;
	SpmIterOrder *type_path;
	SpmIterOrder *type_ignore;

	Node(uint32_t depth_)
	: depth(depth_) {
		uint32_t i;
		this->type_path = new SpmIterOrder[XFORM_MAX];
		this->type_ignore = new SpmIterOrder[XFORM_MAX];
		for (i=0; i<((uint32_t) XFORM_MAX); i++) {
			this->type_path[i] = NONE;
			this->type_ignore[i] = NONE;
		}
	}
	~Node() {}
	void Insert(SpmIterOrder type);
	void Ignore(SpmIterOrder type);
	Node MakeChild(SpmIterOrder type);
	void PrintNode();
};

} // end of csx namespace

#endif /* CSX_DRLE_H__ */
