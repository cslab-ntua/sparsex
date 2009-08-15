#include "spm.h"
#include "mmf.h"

#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <boost/function.hpp>
#include <boost/foreach.hpp>
#define FOREACH BOOST_FOREACH
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include "../../debug/debug_user.h"

extern "C" {
	#include "dynarray.h"
}

namespace bll = boost::lambda;
using namespace csx;

//
// SPM::Builder hide the detaisl of dynarray
//

SPM::Builder::Builder(SPM *_spm, uint64_t elems_nr, uint64_t rows_nr): spm(_spm)
{
	uint64_t *rowptr;

	this->da_elems = dynarray_create(sizeof(SpmRowElem), elems_nr ? elems_nr : 512);
	this->da_rowptr = dynarray_create(sizeof(uint64_t), rows_nr ? rows_nr : 512);

	rowptr = (uint64_t *)dynarray_alloc(this->da_rowptr);
	*rowptr = 0;
}

SPM::Builder::~Builder()
{
	// check if finalized has been called
	assert(this->da_elems == NULL);
	assert(this->da_rowptr == NULL);
}

SpmRowElem *SPM::Builder::AllocElem()
{
	return (SpmRowElem *)dynarray_alloc(this->da_elems);
}

SpmRowElem *SPM::Builder::AllocElems(uint64_t nr)
{
	return (SpmRowElem *)dynarray_alloc_nr(this->da_elems, nr);
}

uint64_t SPM::Builder::getElemsCnt()
{
	return dynarray_size(this->da_elems);
}

void SPM::Builder::newRow(uint64_t rdiff)
{
	uint64_t *rowptr;
	uint64_t elems_cnt;

	elems_cnt = getElemsCnt();
	rowptr = (uint64_t *)dynarray_alloc_nr(this->da_rowptr, rdiff);

	for (uint64_t i=0; i < rdiff; i++)
		rowptr[i] = elems_cnt;
}

void SPM::Builder::Finalize()
{
	// Add marker for last row, if needed
	uint64_t *last_rowptr;
	last_rowptr = (uint64_t *)dynarray_get_last(this->da_rowptr);
	if (*last_rowptr != dynarray_size(this->da_elems))
		newRow();

	// release old data from spm, if needed
	if (spm->elems__)
		free(spm->elems__);
	if (spm->rowptr__)
		free(spm->rowptr__);

	// fill new data in spm
	spm->elems_size__ = dynarray_size(da_elems);
	spm->elems__ = (SpmRowElem *)dynarray_destroy(da_elems);
	spm->rowptr_size__ = dynarray_size(da_rowptr);
	spm->rowptr__ = (uint64_t *)dynarray_destroy(da_rowptr);

	// done!
	this->da_elems = this->da_rowptr = NULL;
}

namespace csx {

//
// Set Elems: build spm via point iterators
//

void mk_row_elem(const CooElem &p, SpmRowElem *ret)
{
	ret->x = p.x;
	ret->val = p.val;
	ret->pattern = NULL;
}

void mk_row_elem(const SpmCooElem &p, SpmRowElem *ret)
{
	ret->x = p.x;
	ret->val = p.val;
	ret->pattern = (p.pattern == NULL) ? NULL : (p.pattern)->clone();
}

void mk_row_elem(const SpmRowElem &p, SpmRowElem *ret)
{
	ret->x = p.x;
	ret->val = p.val;
	ret->pattern = (p.pattern == NULL) ? NULL : (p.pattern)->clone();
}

} // end namespace csx

template <typename IterT>
uint64_t SPM::SetElems(IterT &pi, const IterT &pnts_end,
                       const unsigned long limit,
                       const uint64_t elems_nr, const uint64_t rows_nr)
{
	SpmRowElem *elem;
	uint64_t row_prev, row;
	SPM::Builder *SpmBld;

	SpmBld = new SPM::Builder(this, elems_nr, rows_nr);

	row_prev = 1;
	for (; pi != pnts_end; ++pi){
		row = (*pi).y;
		if (row != row_prev){
			assert(row > row_prev);
			if (limit && SpmBld->getElemsCnt() >= limit)
				break;
			SpmBld->newRow(row - row_prev);
			row_prev = row;
		}
		// add new element
		elem = SpmBld->AllocElem();
		mk_row_elem(*pi, elem);
	}

	SpmBld->Finalize();
	delete SpmBld;

	return this->elems_size__;
}

SPM *SPM::loadMMF_mt(std::istream &in, const long nr)
{
	SPM *ret, *spm;
	MMF mmf(in);
	MMF::iterator iter = mmf.begin();
	MMF::iterator iter_end = mmf.end();
	long limit, cnt, row_start;

	ret = new SPM[nr];

	row_start = limit = cnt = 0;
	for (long i=0; i<nr; i++){
		spm = ret + i;
		limit = (mmf.nnz - cnt) / (nr - i);
		spm->nnz = spm->SetElems(iter, iter_end, limit);
		spm->nrows = spm->rowptr_size__ - 1;
		spm->ncols = mmf.ncols;
		spm->row_start = row_start;
		row_start += spm->nrows;
		spm->type = HORIZONTAL;
		cnt += spm->nnz;
	}
	assert((uint64_t)cnt == mmf.nnz);

	return ret;
}


SPM *SPM::loadMMF_mt(const char *mmf_file, const long nr)
{
	SPM *ret;
	std::ifstream mmf;

	mmf.open(mmf_file);
	ret = loadMMF_mt(mmf, nr);
	mmf.close();

	return ret;
}

SPM *SPM::loadMMF(std::istream &in)
{
	return loadMMF_mt(in, 1);
}

SPM *SPM::loadMMF(const char *mmf_file)
{
	SPM *ret;

	std::ifstream mmf;
	mmf.open(mmf_file);
	ret = loadMMF(mmf);
	mmf.close();

	return ret;
}

namespace csx {


std::ostream &operator<<(std::ostream &os, const Pattern::StatsVal &stats)
{
	os << "nnz: " << stats.nnz;
	return os;
}

std::ostream &operator<<(std::ostream &os, const Pattern &p)
{
	os << " (";
	p.print_on(os);
	os << " type:" << p.type << ") ";
	return os;
}

std::ostream &operator<<(std::ostream &out, CooElem p)
{
	out << "(" << std::setw(2) << p.y << "," << std::setw(2) << p.x << ")";
	return out;
}

std::ostream &operator<<(std::ostream &out, const SpmCooElem e)
{
	out << static_cast<CooElem>(e);
	if (e.pattern != NULL){
		out << "->[" << *(e.pattern) << "]";
		out << " vals:{ ";
		for (int i=0; i<e.pattern->getSize(); i++){
			out << e.vals[i] << " ";
		}
		out << "}";
	} else {
		out << "v: " << e.val;
	}
	return out;
}

std::ostream &operator<<(std::ostream &out, const SpmRowElem &elem)
{
	out << "x:" << elem.x;
	if (elem.pattern){
		out << *(elem.pattern);
		out << " vals:{ ";
		for (int i=0; i<elem.pattern->getSize(); i++){
			out << elem.vals[i] << " ";
		}
		out << "}";
	} else {
		out << " v:" << elem.val;
	}
	return out;
}

std::ostream &operator<<(std::ostream &out, SPM::PntIter pi)
{
	out << "<" << std::setw(2) << pi.row_idx << "," << std::setw(2) << pi.elm_idx << ">";
	return out;
}

// mappings for vertical transformation
static inline void pnt_map_V(CooElem &src, CooElem &dst)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	dst.x = src_y;
	dst.y = src_x;
}
static inline void pnt_rmap_V(CooElem &src, CooElem &dst)
{
	pnt_map_V(src, dst);
}

static inline void pnt_map_D(const CooElem &src, CooElem &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	assert(nrows + src_x - src_y > 0);
	dst.y = nrows + src_x - src_y;
	dst.x = (src_x < src_y) ? src_x : src_y;
}
static inline void pnt_rmap_D(const CooElem &src, CooElem &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	if (src_y < nrows) {
		dst.x = src_x;
		dst.y = nrows + src_x - src_y;
	} else {
		dst.y = src_x;
		dst.x = src_y + src_x - nrows;
	}
}

static inline void pnt_map_rD(const CooElem &src, CooElem &dst, uint64_t ncols)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	uint64_t dst_y;
	dst.y = dst_y = src_x + src_y - 1;
	dst.x = (dst_y <= ncols) ? src_y : ncols + 1 - src_x;
}

static inline void pnt_rmap_rD(const CooElem &src, CooElem &dst, uint64_t ncols)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	if (src_y <= ncols){
		dst.y = src_x;
		dst.x = src_y + 1 - src_x;
	} else {
		dst.y = src_x + src_y - ncols;
		dst.x = ncols + 1 - src_x;
	}
	dst.y = src_y - dst.x + 1;
}

} // end of csx namespace

// Spm Row Elements iteration
SpmRowElem *SPM::rbegin(uint64_t ridx)
{
	assert(ridx < this->rowptr_size__ - 1);
	return &this->elems__[this->rowptr__[ridx]];
}

SpmRowElem *SPM::rend(uint64_t ridx)
{
	assert(ridx < this->rowptr_size__ - 1);
	return &this->elems__[this->rowptr__[ridx + 1]];
}


//
// Spm Point Iteration
//

SPM::PntIter::PntIter(): spm(NULL), row_idx(0), elm_idx(0) { }

SPM::PntIter::PntIter(SPM *_spm, uint64_t ridx) : spm(_spm), row_idx(ridx)
{
	uint64_t *rp = this->spm->rowptr__;
	uint64_t rp_size = this->spm->rowptr_size__;

	assert(ridx < rp_size);
	while (ridx < rp_size && rp[ridx] == rp[ridx+1])
		ridx++;
	this->row_idx = ridx;
	this->elm_idx = rp[ridx];
}

bool SPM::PntIter::operator==(const PntIter &pi)
{
	return (spm = pi.spm)
	       && (row_idx == pi.row_idx)
		   && (elm_idx == pi.elm_idx);
}

bool SPM::PntIter::operator!=(const PntIter &pi)
{
	return !(*this == pi);
}

void SPM::PntIter::operator++()
{
	uint64_t *rp = this->spm->rowptr__;
	uint64_t rp_size = this->spm->rowptr_size__;

	assert(this->elm_idx < this->spm->elems_size__);
	assert(this->row_idx < rp_size);

	this->elm_idx++;
	while (this->row_idx < rp_size && rp[this->row_idx +1] == this->elm_idx)
		this->row_idx++;
}

SpmCooElem SPM::PntIter::operator*()
{
	SpmCooElem ret;
	SpmRowElem *e;
	Pattern *p;

	ret.y = this->row_idx + 1;
	e = this->spm->elems__ + this->elm_idx;
	ret.x = e->x;
	ret.val = e->val;
	p = e-> pattern;
	ret.pattern = (p == NULL) ? NULL : p->clone();
	return ret;
}

SPM::PntIter SPM::points_begin(uint64_t ridx)
{
	return PntIter(this, ridx);
}

SPM::PntIter SPM::points_end(uint64_t ridx)
{
	if (ridx == 0){
		ridx = this->getNrRows();
	}
	return PntIter(this, ridx);
}

static inline bool elem_cmp_less(const SpmCooElem &e0,
                                 const SpmCooElem &e1)
{
	int ret;
	ret = CooCmp(static_cast<CooElem>(e0), static_cast<CooElem>(e1));
	return (ret < 0);
}

inline TransformFn SPM::getRevXformFn(SpmIterOrder type)
{
	boost::function<void (CooElem &p)> ret;
	switch(type) {
		case HORIZONTAL:
		break;

		case VERTICAL:
		ret = bll::bind(pnt_rmap_V, bll::_1, bll::_1);
		break;

		case DIAGONAL:
		ret = bll::bind(pnt_rmap_D, bll::_1, bll::_1, this->nrows);
		break;

		case REV_DIAGONAL:
		ret = bll::bind(pnt_rmap_rD, bll::_1, bll::_1, this->ncols);
		break;

		default:
		std::cerr << "Unknown type: " << type << std::endl;
		assert(false);
	}
	return ret;
}

inline TransformFn SPM::getXformFn(SpmIterOrder type)
{
	boost::function<void (CooElem &p)> ret;
	switch(type) {
		case VERTICAL:
		ret = bll::bind(pnt_map_V, bll::_1, bll::_1);
		break;

		case DIAGONAL:
		ret = bll::bind(pnt_map_D, bll::_1, bll::_1, this->nrows);
		break;

		case REV_DIAGONAL:
		ret = bll::bind(pnt_map_rD, bll::_1, bll::_1, this->ncols);
		break;

		case HORIZONTAL:
		ret = NULL;
		break;

		default:
		assert(false);
	}
	return ret;
}

inline TransformFn SPM::getTransformFn(SpmIterOrder from_type, SpmIterOrder to_type)
{
	boost::function<void (CooElem &p)> xform_fn, rxform_fn;

	rxform_fn = getRevXformFn(from_type);
	xform_fn = getXformFn(to_type);
	if (xform_fn == NULL){
		// to_type is the default type. Just use the rxform_fn
		return rxform_fn;
	}
	if (rxform_fn != NULL){
		// do a double xform: this->type -> HORIZONTAL -> t
		xform_fn = bll::bind(xform_fn, (bll::bind(rxform_fn, bll::_1), bll::_1));
	}
	return xform_fn;
}

void SPM::Transform(SpmIterOrder t, uint64_t rs, uint64_t re)
{
	PntIter p0, pe, p;
	std::vector<SpmCooElem> elems;
	boost::function<void (CooElem &p)> xform_fn;

	if (this->type == t)
		return;

	// Get the appropriate transformation function
	xform_fn = this->getTransformFn(this->type, t);

	// Create a vector of points
	//  In the first version of this we used special iterators that
	//  removed elements, for minimal memory usage.
	//  For now we keep it simple.
	p0 = points_begin(rs);
	pe = points_end(re);
	elems.reserve(this->elems_size__);
	for(p=p0; p != pe; ++p){
		SpmCooElem p_new = SpmCooElem(*p);
		xform_fn(p_new);
		elems.push_back(p_new);
	}

	// Note that if we need minimal memory usage, this part
	// has to change too.
	sort(elems.begin(), elems.end(), elem_cmp_less);
	std::vector<SpmCooElem>::iterator e0, ee;
	e0 = elems.begin();
	ee = elems.end();
	SetElems(e0, ee);
	elems.clear();

	this->type = t;
}



#if 0
void SPM::PrintRows(std::ostream &out)
{
	uint64_t prev_y=1;
	out << prev_y << ": ";
	FOREACH(SpmCooElem &e, this->elems){
		uint64_t y = e.y;
		uint64_t x = e.x;
		if (y != prev_y){
			out << std::endl << y << ": ";
			prev_y = y;
		}
		out << x << " ";
	}
	out << std::endl;
}



void SPM::genStats()
{
}
void PrintTransform(uint64_t y, uint64_t x, TransformFn xform_fn, std::ostream &out)
{
	uint64_t i,j;
	point_t p;
	for (i=1; i <= y; i++){
		for (j=1; j <= x;  j++){
			p.y = i;
			p.x = j;
			xform_fn(p);
			out << p << " ";
		}
		out << std::endl;
	}
}

void PrintDiagTransform(uint64_t y, uint64_t x, std::ostream &out)
{
	boost::function<void (point_t &p)> xform_fn;
	xform_fn = boost::bind(pnt_map_D, _1, _1, y);
	PrintTransform(y, x, xform_fn, out);
}

void PrintRDiagTransform(uint64_t y, uint64_t x, std::ostream &out)
{
	boost::function<void (point_t &p)> xform_fn;
	xform_fn = boost::bind(pnt_map_rD, _1, _1, y);
	PrintTransform(y, x, xform_fn, out);
}

void printXforms()
{
	std::cout << "Diagonal" << std::endl;
	PrintDiagTransform(10, 10, std::cout);
	std::cout << std::endl;
	PrintDiagTransform(5, 10, std::cout);
	std::cout << std::endl;
	PrintDiagTransform(10, 5, std::cout);
	std::cout << std::endl;
	std::cout << "Reverse Diagonal" << std::endl;
	PrintRDiagTransform(10, 10, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(5, 10, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(10, 5, std::cout);
	std::cout << std::endl;
}

void TestTransform(uint64_t y, uint64_t x, TransformFn xform_fn, TransformFn rxform_fn)
{
	uint64_t i,j;
	point_t p0, p1;
	for (i=1; i <= y; i++){
		for (j=1; j <= x;  j++){
			p0.y = i;
			p0.x = j;
			xform_fn(p0);
			p1.y = p0.y;
			p1.x = p0.x;
			rxform_fn(p1);
			if ( (p1.y != i) || (p1.x != j) ){
				std::cout << "Error for " << i << "," << j << std::endl;
				std::cout << "Transformed: " << p0.y << "," << p0.x << std::endl;
				std::cout << "rTransformed:" << p1.y << "," << p1.x << std::endl;
				exit(1);
			}
		}
	}
}

void TestDiagTransform(uint64_t y, uint64_t x)
{
	boost::function<void (point_t &p)> xform_fn, rxform_fn;
	xform_fn = boost::bind(pnt_map_D, _1, _1, y);
	rxform_fn = boost::bind(pnt_rmap_D, _1, _1, y);
	TestTransform(y, x, xform_fn, rxform_fn);
}

void TestRDiagTransform(uint64_t y, uint64_t x)
{
	boost::function<void (point_t &p)> xform_fn, rxform_fn;
	xform_fn = boost::bind(pnt_map_rD, _1, _1, y);
	rxform_fn = boost::bind(pnt_rmap_rD, _1, _1, y);
	TestTransform(y, x, xform_fn, rxform_fn);
}

void TestXforms()
{
	TestDiagTransform(10, 10);
	TestDiagTransform(5, 10);
	TestDiagTransform(10, 5);
	TestRDiagTransform(10, 10);
	TestRDiagTransform(5, 10);
	TestRDiagTransform(10, 5);
}
#endif

void SPM::Print(std::ostream &out)
{
	SPM::PntIter p, p_start, p_end;
	p_start = this->points_begin();
	p_end = this->points_end();
	for (p = p_start; p != p_end; ++p){
		out << " " << (*p);
	}
	out << std::endl;
}

#if 0
void test_drle()
{
	std::vector<int> v_in, delta;
	std::vector< RLE<int> > drle;

	v_in.resize(16);
	for (int i=0; i<16; i++){
		v_in[i] = i;
	}
	FOREACH(int &v, v_in){
		std::cout << v << " ";
	}
	std::cout << std::endl;

	delta = DeltaEncode(v_in);
	FOREACH(int &v, delta){
		std::cout << v << " ";
	}
	std::cout << std::endl;

	drle = RLEncode(delta);
	FOREACH(RLE<int> &v, drle){
		std::cout << "(" << v.val << "," << v.freq << ") ";
	}
	std::cout << std::endl;
}
#endif

#if 0
int main(int argc, char **argv)
{
	SPM obj;
	obj.loadMMF();
	//obj.DRLEncode();
	//obj.Print();
	//std::cout << obj.rows;
	//obj.Print();
	//obj.Draw("1.png");
	//std::cout << "\n";
	//obj.Transform(VERTICAL);
	//std::cout << obj.rows;
	//obj.Print();
	//obj.DRLEncode();
	//std::cout << obj.rows;
}
#endif
