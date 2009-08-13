#include <algorithm>

#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>
#define FOREACH BOOST_FOREACH

#include "spm.h"
#include "drle.h"

namespace bll = boost::lambda;

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
		if (rle.val != curr){
			output.push_back(rle);
			rle.freq = 1;
			rle.val = curr;
		} else {
			rle.freq++;
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
	SPM *Spm;
	std::vector<uint64_t> xs;
	DeltaRLE::Stats stats;

	Spm = this->spm;

	for (uint64_t i=0; i < Spm->getNrRows(); i++){
		for (const SpmRowElem *elem = Spm->rbegin(i); elem != Spm->rend(i); elem++){
			if (elem->pattern == NULL){
				xs.push_back(elem->x);
				continue;
			}
			this->updateStats(xs, stats);
		}
		this->updateStats(xs, stats);
	}

	return stats;
}

// Use to encode a part of a row
//  xs: x values to encode
//  vs: numerical values for the elements
//  newrow: vector to append the encoded elements
void DRLE_Manager::doEncode(std::vector<uint64_t> &xs,
                            std::vector<double> &vs,
                            std::vector<SpmRowElem> &newrow)
{
	uint64_t col; // keep track of the current column
	std::vector< RLE<uint64_t> > rles; // rle elements
	std::set<uint64_t> *deltas_set; // delta values to encode for
	std::vector<double>::iterator vi = vs.begin(); // value iterator
	SpmRowElem elem; // temp element to perform insertions

	// do a delta run-length encoding of the x values
	rles = RLEncode(DeltaEncode(xs));

	// Not all delta rles are to be encoded, only those
	// that are in the ->DeltasToEncode set
	deltas_set = &this->DeltasToEncode[this->spm->type];

	col = 0; // initialize column
	elem.pattern = NULL; // Default inserter (for push_back copies)
	FOREACH(RLE<uint64_t> rle, rles){
		//std::cout << "freq:" << rle.freq << " val:" << rle.val << "\n";
		// create patterns
		if ( deltas_set->find(rle.val) != deltas_set->end() ){
			while (rle.freq >= this->min_limit){
				uint64_t freq;
				SpmRowElem *last_elem;

				freq = std::min(this->max_limit, rle.freq);
				col += rle.val;
				elem.x = col;
				newrow.push_back(elem);

				// get a reference to the last element (avoid unnecessary copies)
				last_elem = &newrow.back();
				// set pattern
				last_elem->pattern = new DeltaRLE(freq, rle.val, this->spm->type);
				// set values
				last_elem->vals = new double[freq];
				std::copy(vi, vi + freq, last_elem->vals);
				vi += freq;

				col += rle.val*(freq - 1);
				rle.freq -= freq;
			}
		}

		// add individual elements
		for (int i=0; i < rle.freq; i++){
			col += rle.val;
			elem.x = col;
			elem.val = *vi++;
			newrow.push_back(elem);
		}
	}

	assert(vi == vs.end());
	xs.clear();
	vs.clear();
}

void DRLE_Manager::EncodeRow(const SpmRowElem *rstart, const SpmRowElem *rend,
                             std::vector<SpmRowElem> &newrow)
{
	std::vector<uint64_t> xs;
	std::vector<double> vs;

	// gather x values into xs vector until a pattern is found
	// and encode them using doEncode()
	for (const SpmRowElem *e = rstart; e < rend; e++){
		if (e->pattern == NULL){
			xs.push_back(e->x);
			vs.push_back(e->val);
			continue;
		}
		if (xs.size() != 0){
			doEncode(xs, vs, newrow);
		}
		newrow.push_back(*e);
	}

	// Encode any remaining elements
	if (xs.size() != 0){
		doEncode(xs, vs, newrow);
	}
}

void DRLE_Manager::Encode(SpmIterOrder type)
{
	SPM *Spm;
	SPM::Builder *SpmBld;
	SpmIterOrder oldtype;
	std::vector<SpmRowElem> new_row;
	uint64_t nr_size;
	SpmRowElem *elems;

	Spm = this->spm;
	if (type == NONE && ((type = this->chooseType()) == NONE) ){
		return;
	}

	// Transform matrix to the desired iteration order
	oldtype = Spm->type;
	Spm->Transform(type);

	// Do the encoding
	SpmBld = new SPM::Builder(Spm);
	for (uint64_t i=0; i < Spm->getNrRows(); i++){

		EncodeRow(Spm->rbegin(i), Spm->rend(i), new_row);

		nr_size = new_row.size();
		if (nr_size > 0){
			elems = SpmBld->AllocElems(nr_size);
			for (uint64_t i=0; i < nr_size; i++){
				mk_row_elem(new_row[i], elems + i);
			}
		}
		new_row.clear();
		SpmBld->newRow();
	}
	SpmBld->Finalize();
	delete SpmBld;

	// Transform matrix to the original iteration order
	Spm->Transform(oldtype);

	this->addIgnore(type);
}

void DRLE_Manager::EncodeAll()
{
	SpmIterOrder type;

	for (;;){
		this->genAllStats();
		this->outStats(std::cerr);
		type = this->chooseType();
		if (type == NONE)
			break;
		std::cerr << "Encode to " << SpmTypesNames[type] << std::endl;
		this->Encode(type);
	}
}

Pattern::Generator *DeltaRLE::generator(CooElem start)
{
	DeltaRLE::Generator *g;
	g = new DeltaRLE::Generator(start, this);
	return g;
}

namespace csx {
void DRLE_OutStats(DeltaRLE::Stats &stats, SPM &spm, std::ostream &os)
{
	DeltaRLE::Stats::iterator iter;
	for (iter=stats.begin(); iter != stats.end(); ++iter){
		os << "    " << iter->first << "-> "
		   << "np:" << iter->second.npatterns
		   << " nnz: " <<  100*((double)iter->second.nnz/(double)spm.nnz) << "%"
		   << " (" << iter->second.nnz << ")";
	}
}
} // end csx namespace

void DRLE_Manager::addIgnore(SpmIterOrder type)
{
	this->xforms_ignore.set(type);
}

void DRLE_Manager::genAllStats()
{
	DeltaRLE::Stats::iterator iter, tmp;
	DeltaRLE::Stats *sp;

	this->stats.clear();
	for (int t=HORIZONTAL; t != XFORM_MAX; t++){
		if (this->xforms_ignore[t])
			continue;

		SpmIterOrder type = SpmTypes[t];
		this->spm->Transform(type);
		this->stats[type] = this->generateStats();
		this->spm->Transform(HORIZONTAL);

		// ** Filter stats
		// From http://www.sgi.com/tech/stl/Map.html:
		// Map has the important property that inserting a new element into a
		// map does not invalidate iterators that point to existing elements.
		// Erasing an element from a map also does not invalidate any
		// iterators, except, of course, for iterators that actually point to
		// the element that is being erased.
		sp = &this->stats[type];
		for (iter = sp->begin(); iter != sp->end(); ){
			tmp = iter++;
			double p = (double)tmp->second.nnz/(double)spm->nnz;
			if (p < this->min_perc){
				sp->erase(tmp);
			} else {
				this->DeltasToEncode[type].insert(tmp->first);
			}
		}
	}
}

// get the number of non-zero elements that can be encoded
// using drle for a specific iteration order (matrix type)
uint64_t DRLE_Manager::getTypeNNZ(SpmIterOrder type)
{
	DeltaRLE::Stats *sp;
	DeltaRLE::Stats::iterator iter;
	uint64_t ret;

	ret = 0;
	if (this->stats.find(type) == this->stats.end())
		return ret;

	sp = &this->stats[type];
	for (iter=sp->begin(); iter != sp->end(); ++iter){
		ret += iter->second.nnz;
	}
	return ret;
}

// choose a type to encode the matrix, based on the stats
// (whichever maximizes getTypeNNZ())
SpmIterOrder DRLE_Manager::chooseType()
{
	SpmIterOrder ret;
	uint64_t max_out;
	DRLE_Manager::StatsMap::iterator iter;

	ret = NONE;
	max_out = 0;
	for (iter=this->stats.begin(); iter != this->stats.end(); ++iter){
		uint64_t out = this->getTypeNNZ(iter->first);
		if (out > max_out){
			max_out = out;
			ret = iter->first;
		}
	}

	return ret;
}

void DRLE_Manager::outStats(std::ostream &os)
{
	DRLE_Manager::StatsMap::iterator iter;
	for (iter = this->stats.begin(); iter != this->stats.end(); ++iter){
		os << SpmTypesNames[iter->first] << "\t";
		DRLE_OutStats(iter->second, *(this->spm), os);
		os << std::endl;
	}
}
