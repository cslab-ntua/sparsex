#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <ios>
#include <exception>
#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>
#define FOREACH BOOST_FOREACH

#include "spm.h"
#include "drle.h"

namespace bll = boost::lambda;

using namespace csx;

bool debug = false;

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
std::vector<RLE<typename T::iterator::value_type> >
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

DRLE_Manager::DRLE_Manager(SPM *_spm,
                           long min_limit_,
                           long max_limit_,
                           double min_perc_,
                           uint64_t sort_window_size_,
                           split_alg_t split_type_,
                           double probability)
	: spm(_spm), min_limit(min_limit_), max_limit(max_limit_),
      min_perc(min_perc_), sort_window_size(sort_window_size_),
      split_type(split_type_) {
    // These are delimiters, ignore them by default.
    addIgnore(BLOCK_TYPE_START);
    addIgnore(BLOCK_COL_START);
    addIgnore(BLOCK_TYPE_END);
    addIgnore(BLOCK_ROW_TYPE_NAME(1));
    addIgnore(BLOCK_COL_TYPE_NAME(1));

    check_and_set_sorting();
    if (sort_windows) {
        compute_sort_splits();

        // Initialize sampling stuff
        srand48(0);
        set_sampling_probability(probability);
    }
}


void DRLE_Manager::check_and_set_sorting()
{
    switch (split_type) {
    case DRLE_Manager::SPLIT_BY_ROWS:
        do_check_sort_by_rows();
        break;
    case DRLE_Manager::SPLIT_BY_NNZ:
        do_check_sort_by_nnz();
        break;
    default:
        throw new std::invalid_argument("Unknown split algorithm");
    }
}

void DRLE_Manager::do_check_sort_by_rows()
{
    if (sort_window_size > spm->getNrRows())
        throw new std::invalid_argument("Invalid sort window");

    if (sort_window_size == 0 ||
        sort_window_size == spm->getNrRows())
        sort_windows = false;
    else
        sort_windows = true;
}

void DRLE_Manager::do_check_sort_by_nnz()
{
    if (sort_window_size > spm->elems_size__)
        throw new std::invalid_argument("Invalid sort window");

    if (sort_window_size == 0 ||
        sort_window_size == spm->elems_size__)
        sort_windows = false;
    else
        sort_windows = true;
}

void DRLE_Manager::compute_sort_splits()
{
    switch (split_type) {
    case DRLE_Manager::SPLIT_BY_ROWS:
        do_compute_sort_splits_by_rows();
        break;
    case DRLE_Manager::SPLIT_BY_NNZ:
        do_compute_sort_splits_by_nnz();
        break;
    default:
        throw new std::invalid_argument("Unknown split algorithm");
    }
}

void DRLE_Manager::do_compute_sort_splits_by_rows()
{
    uint64_t nr_rows = spm->getNrRows();
    uint64_t i;
    for (i = 0; i <= nr_rows; i += sort_window_size)
        sort_splits.push_back(i);
    if (i > nr_rows && i - nr_rows < sort_window_size / 2)
        sort_splits.push_back(nr_rows);
    else {
        // the last window is too short; merge it with the previous
        sort_splits.pop_back();
        sort_splits.push_back(nr_rows);
    }
}

void DRLE_Manager::do_compute_sort_splits_by_nnz()
{
    uint64_t nzeros_cnt;
    uint64_t nr_rows = spm->getNrRows();

    nzeros_cnt = 0;
    sort_splits.push_back(0);
    for (uint64_t i = 0; i < nr_rows; ++i) {
        uint64_t new_nzeros_cnt =
            nzeros_cnt + spm->rowptr__[i+1] - spm->rowptr__[i];
        if (new_nzeros_cnt < sort_window_size) {
            nzeros_cnt = new_nzeros_cnt;
        } else {
            // new split
            sort_splits.push_back(i+1);
            nzeros_cnt = 0;
        }
    }

    if (sort_splits.back() < nr_rows)
        sort_splits.push_back(nr_rows);
}

//
// Iterate over stats of `type' and multiply them by `factor'
// 
void DRLE_Manager::correct_stats(SpmIterOrder type, double factor)
{
    if (stats.find(type) == stats.end())
        return;

//    std::cout << "correcting factor: " << factor << std::endl;
    DeltaRLE::Stats *l_stats = &stats[type];
    DeltaRLE::Stats::iterator   iter;
    for (iter = l_stats->begin(); iter != l_stats->end(); ++iter) {
        iter->second.nnz = (uint64_t) (iter->second.nnz*factor);
        iter->second.npatterns = (long) (iter->second.npatterns*factor);
        assert(iter->second.nnz <= spm->nnz);
    }
}


void DRLE_Manager::updateStats(SpmIterOrder type, DeltaRLE::Stats stats)
{
    if (this->stats.find(type) == this->stats.end()) {
        // first time for this type
        this->stats[type] = stats;
    } else {
        // update stats
        for (DeltaRLE::Stats::const_iterator it = stats.begin();
             it != stats.end(); ++it) {
            this->stats[type][it->first].update(it->second);
        }
    }
}

void DRLE_Manager::updateStats(std::vector<uint64_t> &xs,
                               DeltaRLE::Stats &stats)
{
    updateStats(this->spm, xs, stats);
}

void DRLE_Manager::updateStats(SPM *spm, std::vector<uint64_t> &xs,
                               DeltaRLE::Stats &stats)
{
	std::vector< RLE<uint64_t> > rles;
    uint64_t    block_align = isBlockType(spm->type);
    if (block_align) {
        DRLE_Manager::updateStatsBlock(xs, stats, block_align);
        return;
    }

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

void DRLE_Manager::updateStatsBlock(std::vector<uint64_t> &xs,
                                    DeltaRLE::Stats &stats,
                                    uint64_t block_align)
{
	std::vector< RLE<uint64_t> > rles;
	
	assert(block_align);
	if (xs.size() == 0)
		return;

	rles = RLEncode(DeltaEncode(xs));
    	uint64_t unit_start = 0;
        FOREACH(RLE<uint64_t> &rle, rles){
        	unit_start += rle.val;
            // printf("(v,f,u) = (%ld,%ld,%ld)\n", rle.val, rle.freq, unit_start);
        	if (rle.val == 1) {
        		// Start of the real block is at `unit_start - 1' with
                // one-based indexing. When computing the `%' we need
                // zero-based indexing.
                uint64_t nr_elem = rle.freq + 1;
                uint64_t skip_front = (unit_start == 1) ? 0 :
                    (unit_start - 2) % block_align;
                if (nr_elem > skip_front)
                    nr_elem -= skip_front;
                else
                    nr_elem = 0;
                uint64_t other_dim = nr_elem / (uint64_t) block_align;
                if (other_dim >= 2) {
                    stats[other_dim].nnz += other_dim * block_align;
                    stats[other_dim].npatterns++;
                }
        	}
        	unit_start += rle.val*(rle.freq - 1);
        }
        xs.clear();
}

DeltaRLE::Stats DRLE_Manager::generateStats(uint64_t rs, uint64_t re)
{
    return generateStats(this->spm, rs, re);
}

DeltaRLE::Stats DRLE_Manager::generateStats(SPM *Spm, uint64_t rs, uint64_t re)
{
	std::vector<uint64_t> xs;
	DeltaRLE::Stats stats;

	for (uint64_t i=rs; i < re; i++){
		for (const SpmRowElem *elem = Spm->rbegin(i);
             elem != Spm->rend(i); elem++){
			if (elem->pattern == NULL){
				xs.push_back(elem->x);
				continue;
			}
			this->updateStats(Spm, xs, stats);
		}
		this->updateStats(Spm, xs, stats);
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
	uint64_t col; 						// keep track of the current column
	std::vector< RLE<uint64_t> > rles; 			// rle elements
	std::set<uint64_t> *deltas_set; 			// delta values to encode for
	std::vector<double>::iterator vi = vs.begin(); 		// value iterator
	SpmRowElem elem; 					// temp element to perform insertions
    
    	if (isBlockType(this->spm->type)) {
        	doEncodeBlock(xs, vs, newrow);
        	return;
    	}

	// do a delta run-length encoding of the x values
	rles = RLEncode(DeltaEncode(xs));

	// Not all delta rles are to be encoded, only those
	// that are in the ->DeltasToEncode set
	deltas_set = &this->DeltasToEncode[this->spm->type];

	col = 0; // initialize column
	elem.pattern = NULL; // Default inserter (for push_back copies)
	FOREACH(RLE<uint64_t> rle, rles){
		// create patterns
       		//std::cout << "freq:" << rle.freq << " val:" << rle.val << "\n";
		if (deltas_set->find(rle.val) != deltas_set->end()){
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
                last_elem->pattern = new DeltaRLE(freq, rle.val,
                                                  this->spm->type);
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

//static uint64_t nr_lines = 0;
void DRLE_Manager::doEncodeBlock(std::vector<uint64_t> &xs,
                                 std::vector<double> &vs,
                                 std::vector<SpmRowElem> &newrow)
{
	uint64_t col; // keep track of the current column
	std::vector< RLE<uint64_t> > rles; // rle elements
	std::set<uint64_t> *deltas_set; // delta values to encode for
	std::vector<double>::iterator vi = vs.begin(); // value iterator
	SpmRowElem elem; // temp element to perform insertions
    	SpmRowElem *last_elem;
    
	// do a delta run-length encoding of the x values
	rles = RLEncode(DeltaEncode(xs));

	// Not all delta rles are to be encoded, only those
	// that are in the ->DeltasToEncode set
	deltas_set = &this->DeltasToEncode[this->spm->type];

    int block_align = isBlockType(this->spm->type);
    assert(block_align);

	col = 0; // initialize column
	elem.pattern = NULL; // Default inserter (for push_back copies)
	FOREACH(RLE<uint64_t> rle, rles){
	
	// create patterns
        //std::cout << "freq:" << rle.freq << " val:" << rle.val << "\n";

        col += rle.val;
        uint64_t skip_front, skip_back, nr_elem;
        if (col == 1) {
            skip_front = 0;
            nr_elem = rle.freq;
        } else {
            skip_front = (col - 2) % block_align;
            if (skip_front != 0)
                skip_front  = block_align - skip_front;

            nr_elem = rle.freq + 1;
        }
        
        if (nr_elem > skip_front)
            nr_elem -= skip_front;
        else
            nr_elem = 0;

        skip_back = nr_elem % block_align;

        if (nr_elem > skip_back)
            nr_elem -= skip_back;
        else
            nr_elem = 0;

		if (rle.val == 1 &&
            deltas_set->find(nr_elem / block_align) != deltas_set->end() &&
            nr_elem >= (uint64_t) 2*block_align) {

            uint64_t    rle_start;
            // We have a new block RLE
            if (col != 1) {
                rle_start = col - 1;
                // we need to annex the previous element
                newrow.pop_back();
                --vi;
//                std::cout << "popped: " << *vi << std::endl << std::endl;
            } else {
                rle_start = col;
            }

            // Add elements skipped from start
            for (uint64_t i = 0; i < skip_front; i++) {
                elem.x = rle_start + i;
                elem.val = *vi++;
                newrow.push_back(elem);
            }

            // Align max_limit 
            uint64_t max_limit =
                (this->max_limit / (2*block_align))*(2*block_align);

            uint64_t nr_blocks = nr_elem / max_limit;
            uint64_t nr_elem_block = std::min(max_limit, nr_elem);

            if (nr_blocks == 0)
                nr_blocks = 1;
            else
                // Adjust skip_back with leftover items from max_limit alignment
                skip_back += nr_elem - nr_elem_block*nr_blocks;
            
            for (uint64_t i = 0; i < nr_blocks; i++) {
                // Add the blocks
                elem.x = rle_start + skip_front + i*nr_elem_block;
                newrow.push_back(elem);
                last_elem = &newrow.back();
//                if (debug)
//                     std::cout << "nl:" << nr_lines << " rle_start:" << rle_start << " col:" << col << " sf:" << skip_front << " sb:" << skip_back << " nr_elem:" << nr_elem_block << std::endl;

                last_elem->pattern =
                    new BlockRLE(nr_elem_block,
                                 nr_elem_block / block_align, this->spm->type);
                last_elem->vals = new double[nr_elem_block];
                std::copy(vi, vi + nr_elem_block, last_elem->vals);
                vi += nr_elem_block;
            }

            // Add the remaining elements
            for (uint64_t i = 0; i < skip_back; i++) {
                elem.x = rle_start + skip_front + nr_elem_block*nr_blocks + i;
                elem.val = *vi++;
                newrow.push_back(elem);
            }

        } else {
            // add individual elements
            for (int i = 0; i < rle.freq; i++) {
                elem.x = col + i*rle.val;
                elem.val = *vi++;
                newrow.push_back(elem);
            }
        }

        col += rle.val * (rle.freq - 1);
	}

//    std::cerr << vs.end() - vi << std::endl;
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
	for (const SpmRowElem *e = rstart; e < rend; ++e){
		if (e->pattern == NULL){
			xs.push_back(e->x);
			vs.push_back(e->val);
			continue;
		}
		if (xs.size() != 0){
			doEncode(xs, vs, newrow);
		}
		newrow.push_back(*e);
		if (e->pattern)
			delete e->pattern;
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

	// Do the encoding
	SpmBld = new SPM::Builder(Spm);

    // Transform matrix to the desired iteration order
    oldtype = Spm->type;
    Spm->Transform(type);

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
	SpmIterOrder type = NONE;
	StatsMap::iterator iter;

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

void DRLE_Manager::ignoreAll()
{
	this->xforms_ignore.set();
}

void DRLE_Manager::removeIgnore(SpmIterOrder type)
{
	// the following types are always ignored
	if (type <= NONE ||
	    type == BLOCK_TYPE_START ||
	    type == BLOCK_ROW_TYPE_NAME(1) ||
	    type == BLOCK_COL_START ||
	    type == BLOCK_COL_TYPE_NAME(1) ||
	    type == BLOCK_TYPE_END ||
	    type >= XFORM_MAX)
		return;
	this->xforms_ignore.reset(type);
}

void DRLE_Manager::removeAll()
{
    for (int t = NONE; t < XFORM_MAX; t++)
        this->xforms_ignore.reset(t);
}

void DRLE_Manager::genAllStats()
{
	DeltaRLE::Stats::iterator iter, tmp;
	DeltaRLE::Stats *sp;

//     this->addIgnore(HORIZONTAL);
//     this->addIgnore(VERTICAL);
//     this->addIgnore(DIAGONAL);
//     this->addIgnore(REV_DIAGONAL);

//     this->addIgnore(BLOCK_ROW_TYPE_NAME(2));
//     this->addIgnore(BLOCK_ROW_TYPE_NAME(3));
//     this->addIgnore(BLOCK_ROW_TYPE_NAME(4));
//     this->addIgnore(BLOCK_ROW_TYPE_NAME(5));
//     this->addIgnore(BLOCK_ROW_TYPE_NAME(6));
//     this->addIgnore(BLOCK_ROW_TYPE_NAME(7));
//     this->addIgnore(BLOCK_ROW_TYPE_NAME(8));

//     this->addIgnore(BLOCK_COL_TYPE_NAME(2));
//     this->addIgnore(BLOCK_COL_TYPE_NAME(3));
//     this->addIgnore(BLOCK_COL_TYPE_NAME(4));
//     this->addIgnore(BLOCK_COL_TYPE_NAME(5));
//     this->addIgnore(BLOCK_COL_TYPE_NAME(6));
//     this->addIgnore(BLOCK_COL_TYPE_NAME(7));
//     this->addIgnore(BLOCK_COL_TYPE_NAME(8));

	this->stats.clear();

	for (int t=HORIZONTAL; t != XFORM_MAX; t++){
		
		if (this->xforms_ignore[t])
			continue;

		SpmIterOrder type = SpmTypes[t];

        std::cout << "Checking for " << SpmTypesNames[t] << std::endl;
        uint64_t samples_cnt = 0;
		if (sort_windows) {
            sort_split_iterator iter;
            for (iter = sort_splits.begin();
                 iter != sort_splits.end() - 1; ++iter) {
                if (drand48() < 1 - sampling_probability)
                    continue;
                ++samples_cnt;
                uint64_t window_size = *(iter + 1) - *iter;
                SPM *window = this->spm->getWindow(*iter, window_size);
                window->Transform(type);
                DeltaRLE::Stats w_stats = generateStats(window, 0,
                                                        window->getNrRows());
                updateStats(type, w_stats);
                window->Transform(HORIZONTAL);
                this->spm->putWindow(window);
                delete window;
            }

            // std::cout << "nr_splits = " << sort_splits.size() - 1 << ", "
            //           << "nr_samples = " << samples_cnt << std::endl;
            correct_stats(type,
                          (sort_splits.size() - 1) / ((double) samples_cnt));

        } else {
            this->spm->Transform(type);
            this->stats[type] = this->generateStats(0, this->spm->getNrRows());
            this->spm->Transform(HORIZONTAL);
        }

		// ** Filter stats
		// From http://www.sgi.com/tech/stl/Map.html:
		// Map has the important property that inserting a new element into a
		// map does not invalidate iterators that point to existing elements.
		// Erasing an element from a map also does not invalidate any
		// iterators, except, of course, for iterators that actually point to
		// the element that is being erased.
		sp = &this->stats[type];
		for (iter = sp->begin(); iter != sp->end(); ) {
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

//
// Gets a score for each type. This might be used for choosing an encoding.
// 
uint64_t DRLE_Manager::getTypeScore(SpmIterOrder type)
{
	DeltaRLE::Stats *sp;
	DeltaRLE::Stats::iterator iter;
	uint64_t ret;

	ret = 0;
	if (this->stats.find(type) == this->stats.end())
		return ret;
	sp = &this->stats[type];
	uint64_t nr_nzeros_encoded = 0;
	uint64_t nr_patterns = 0;
	for (iter=sp->begin(); iter != sp->end(); ++iter){
		nr_nzeros_encoded += iter->second.nnz;
		nr_patterns += iter->second.npatterns;
	}
	//ret = this->spm->nnz - (nr_patterns + this->spm->nnz - nr_nzeros_encoded);
	ret = nr_nzeros_encoded - nr_patterns;
	return ret;
}

// choose a type to encode the matrix, based on the stats
// (whichever maximizes getTypeScore())
SpmIterOrder DRLE_Manager::chooseType()
{
	SpmIterOrder ret;
	uint64_t max_out;
	DRLE_Manager::StatsMap::iterator iter;

	ret = NONE;
	max_out = 0;
	for (iter=this->stats.begin(); iter != this->stats.end(); ++iter){
		uint64_t out = this->getTypeScore(iter->first);
		if (out == 0) {
			this->addIgnore(iter->first);
		} else if (out > max_out){
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

void DRLE_Manager::doDecode(const SpmRowElem *elem, std::vector<SpmRowElem> &newrow)
{
	long int i;
	uint64_t cur_x;
	SpmRowElem new_elem;	

	new_elem.pattern = NULL;
	cur_x = elem->x;
	for (i=0; i<elem->pattern->getSize(); i++) {
		new_elem.x = cur_x;
		new_elem.val = elem->vals[i];
		newrow.push_back(new_elem);
		cur_x = elem->pattern->getNextX(cur_x);
	}
	delete elem->pattern;
	delete elem->vals;
}

void DRLE_Manager::DecodeRow(const SpmRowElem *rstart, const SpmRowElem *rend, std::vector<SpmRowElem> &newrow)
{
	for (const SpmRowElem *e = rstart; e < rend; e++){
		if (e->pattern != NULL && e->pattern->type == this->spm->type) {
			doDecode(e, newrow);
		}
		else {
			newrow.push_back(*e);
			if (e->pattern)
				delete e->pattern;
		}
	}
}

void DRLE_Manager::Decode(SpmIterOrder type)
{
	SPM *Spm;
	SPM::Builder *SpmBld;
	SpmIterOrder oldtype;
	std::vector<SpmRowElem> new_row;
	uint64_t nr_size;
	SpmRowElem *elems;
	
	Spm = this->spm;
	if (type == NONE) {
		return;
	}

	// Transform matrix to the desired iteration order
	oldtype = Spm->type;
	Spm->Transform(type);

	// Do the decoding
	SpmBld = new SPM::Builder(Spm);
	for (uint64_t i=0; i < Spm->getNrRows(); i++){
		DecodeRow(Spm->rbegin(i), Spm->rend(i), new_row);
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
	this->removeIgnore(type);
}

void Node::Insert(SpmIterOrder type) {
	uint32_t i=0;
	
	while (this->type_path[i] != NONE) {
		assert(this->type_path[i] != type);
		i++;	
	}
	this->type_path[i] = type;
}

void Node::Ignore(SpmIterOrder type) {
	uint32_t i=0;
	
	while (this->type_ignore[i] != NONE) {
		assert(this->type_ignore[i] != type);
		i++;	
	}
	this->type_ignore[i] = type;
}

void Node::PrintNode() {
	uint32_t i=0;
	SpmIterOrder type;
	
	//printf("Depth: %d\n",this->depth);
	while ((type = this->type_path[i]) != NONE) {
		if (i != 0)
			std::cout << ",";
		std::cout << type;
		i++;
	}
	std::cout << std::endl;
}

Node Node::MakeChild(SpmIterOrder type) {
	Node new_node = Node(this->depth+1);

	for (uint32_t i=0; i<((uint32_t) XFORM_MAX); i++) {
		new_node.type_path[i] = this->type_path[i];
		new_node.type_ignore[i] = this->type_ignore[i];
	}
	new_node.Insert(type);
	return new_node;
}

void DRLE_Manager::MakeEncodeTree()
{
	uint32_t i,j;
	uint32_t count=0;
	StatsMap::iterator iter;
	std::vector<Node> nodes;
	std::vector<Node>::iterator it,it2;
	
	nodes.push_back(Node(0));
	while (nodes.size() != 0) {
		i = 0;
		while (nodes[0].type_ignore[i] != NONE) {
			this->addIgnore(nodes[0].type_ignore[i]);
			i++;
		}
		this->genAllStats();
		//this->outStats(std::cout);
		i = 0;
		while (nodes[0].type_ignore[i] != NONE) {
			this->removeIgnore(nodes[0].type_ignore[i]);
			i++;
		}
		i = 0;
		for (iter=this->stats.begin(); iter != this->stats.end(); ++iter) {
			if (this->getTypeScore(iter->first) == 0) {
				nodes[0].Ignore(iter->first);
			}
		}
		it = nodes.begin();
		for (iter=this->stats.begin(); iter != this->stats.end(); ++iter) {
			if (this->getTypeScore(iter->first) > 0) {
				it++;
				it = nodes.insert(it,nodes[0].MakeChild(iter->first));
			}
		}
		if (it != nodes.begin()) {
			//std::cout << "Encode to " << SpmTypesNames[nodes[1].type_path[nodes[1].depth-1]] << std::endl;
			this->Encode(nodes[1].type_path[nodes[1].depth-1]);
		}
		else {
			nodes[0].PrintNode();
			count++;
			if (nodes.size() != 1) {
				i = 0;
				while (nodes[0].type_path[i] == nodes[1].type_path[i])
					i++;
				for (j=1; j<=nodes[0].depth-i; j++) {
					//std::cout << "Decode from " << SpmTypesNames[nodes[0].type_path[nodes[0].depth-j]] << std::endl;
					this->Decode(nodes[0].type_path[nodes[0].depth-j]);
				}
				while (i < nodes[1].depth) {
					//std::cout << "Encode to " << SpmTypesNames[nodes[1].type_path[i]] << std::endl;
					this->Encode(nodes[1].type_path[i]);
					i++;
				}
			}
			else {
				for (j=1; j<=nodes[0].depth; j++) {
					//std::cout << "Decode from " << SpmTypesNames[nodes[0].type_path[nodes[0].depth-j]] << std::endl;
					this->Decode(nodes[0].type_path[nodes[0].depth-j]);
				}
			}
		}
		delete nodes[0].type_path;
		delete nodes[0].type_ignore;
		nodes.erase(nodes.begin());
	}
	std::cout << "Tree has " << count << " possible paths" << std::endl;
}

void DRLE_Manager::EncodeSerial() {
	
	const char *xform_orig = getenv("XFORM_CONF");
	if (xform_orig) {
		int len = strlen(xform_orig) + 1;
		char xform_string[len];
		strncpy(xform_string, xform_orig, len);

		int t = atoi(strtok(xform_string, ","));
		for (uint32_t i=0; i<XFORM_MAX; i++)
			this->addIgnore((SpmIterOrder) i);
		this->removeIgnore(static_cast<SpmIterOrder>(t));
		this->genAllStats();
		this->Encode(static_cast<SpmIterOrder>(t));
		char *token;
		while ( (token = strtok(NULL, ",")) != NULL) {
			t = atoi(token);
			for (uint32_t i=0; i<XFORM_MAX; i++)
				this->addIgnore((SpmIterOrder) i);
			this->removeIgnore(static_cast<SpmIterOrder>(t));
			this->genAllStats();			
			this->Encode(static_cast<SpmIterOrder>(t));
		}
	}
}
