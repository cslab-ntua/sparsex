#include "spm.h"
//#include "drle.h"

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
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>


#include "../../debug/debug_user.h"

namespace bll = boost::lambda;

#define FOREACH BOOST_FOREACH

namespace csx {

std::ostream &operator<<(std::ostream &out, CooElem p)
{
	out << "(" << std::setw(2) << p.y << "," << std::setw(2) << p.x << ")";
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

static inline void pnt_map_rD(const CooElem &src, CooElem &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	uint64_t dst_y;
	dst.y = dst_y = src_x + src_y - 1;
	dst.x = (dst_y <= nrows) ? src_x : src_x + nrows - dst_y;
}

static inline void pnt_rmap_rD(const CooElem &src, CooElem &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	if (src_y < nrows){
		dst.x = src_x;
	} else {
		dst.x = src_x + src_y - nrows;
	}
	dst.y = src_y - dst.x + 1;
}


std::ostream &operator<<(std::ostream &out, SpmIdx::PointIter pi)
{
	out << "<" << std::setw(2) << pi.row_idx << "," << std::setw(2) << pi.elm_idx << ">";
	return out;
}


} // end of csx namespace

using namespace csx;

std::istream &operator>>(std::istream &in, SpmIdx &obj)
{
	obj.loadMMF(in);
	return in;
}


SpmIdx::PointIter SpmIdx::points_begin(uint64_t ridx, uint64_t eidx)
{
	return PointIter(ridx, eidx, this);
}

SpmIdx::PointIter SpmIdx::points_end(uint64_t ridx, uint64_t eidx)
{
	if (ridx == 0){
		ridx = this->rows.size();
	}
	return PointIter(ridx, eidx, this);
}

SpmIdx::PointPoper SpmIdx::points_pop_begin()
{
	return PointPoper(0, 0, this);
}

SpmIdx::PointPoper SpmIdx::points_pop_end()
{
	return PointPoper(this->rows.size(), 0, this);
}

static inline bool elem_cmp_less(const SpmCooElem &e0,
                                 const SpmCooElem &e1)
{
	int ret;
	ret = CooCmp(static_cast<CooElem>(e0), static_cast<CooElem>(e1));
	return (ret < 0);
}


static inline void mk_row_elem(const CooElem &p, SpmRowElem &ret)
{
	ret.x = p.x;
	ret.pattern = NULL;
}

static inline void mk_row_elem(const SpmCooElem &p, SpmRowElem &ret)
{
	ret.x = p.x;
	ret.pattern = (p.pattern == NULL) ? NULL : (p.pattern)->clone();
}

template <typename IterT>
void SpmIdx::SetRows(IterT pnts_start, IterT pnts_end)
{
	SpmRowElems *row_elems;
	uint64_t row_prev;
	//std::cout << "SetRows\n";

	//this->rows.resize(this->nrows);
	this->rows.resize(1);
	row_elems = &this->rows.at(0);
	row_prev = 1;

	// XXX
	IterT point;
	for (point = pnts_start; point < pnts_end; ++point){
		uint64_t row = point->y;
		long size;

		if (row != row_prev){
			assert(row > this->rows.size());
			this->rows.resize(row);
			row_elems = &this->rows.at(row - 1);
			row_prev = row;
		}
		size = row_elems->size();
		row_elems->resize(size+1);
		mk_row_elem(*point, (*row_elems)[size]);
	}
	//std::cout << "SetRows [end]\n";
}

void SpmIdx::loadMMF(std::istream &in)
{
	char buff[1024];
	double val;
	uint64_t cnt;
	int ret;
	CooElem *point;
	SpmPoints points;

	// header
	do {
		in.getline(buff, sizeof(buff));
	} while (buff[0] == '#');
	ret = sscanf(buff, "%lu %lu %lu", &this->nrows, &this->ncols, &this->nnz);
	if (ret != 3){
		std::cerr << "mmf header error: sscanf" << std::endl;
		exit(1);
	}

	// body
	type = HORIZONTAL;
	points.resize(this->nnz);
	point = &points[0];
	for (cnt=0; cnt<this->nnz; cnt++){
		in.getline(buff, sizeof(buff));
		ret = sscanf(buff, "%lu %lu %lf", &point->y, &point->x, &val);
		assert(ret == 3);
		point++;
	}
	assert((uint64_t)(point - &points[0]) == this->nnz);
	assert(in.getline(buff, sizeof(buff)).eof());

	SetRows(points.begin(), points.end());
	points.clear();
}

void SpmIdx::loadMMF(const char *mmf_file)
{
	std::ifstream mmf;
	mmf.open(mmf_file);
	this->loadMMF(mmf);
	mmf.close();

}

inline TransformFn SpmIdx::getRevXformFn(SpmIterOrder type)
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
		ret = bll::bind(pnt_rmap_rD, bll::_1, bll::_1, this->nrows);
		break;

		default:
		assert(false);
	}
	return ret;
}

inline TransformFn SpmIdx::getXformFn(SpmIterOrder type)
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
		ret = bll::bind(pnt_map_rD, bll::_1, bll::_1, this->nrows);
		break;

		case HORIZONTAL:
		ret = NULL;

		default:
		assert(false);
	}
	return ret;
}

inline TransformFn SpmIdx::getTransformFn(SpmIterOrder from_type, SpmIterOrder to_type)
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

void SpmIdx::Transform(SpmIterOrder t)
{
	PointPoper p0, pe, p;
	SpmCooElems elems;
	boost::function<void (CooElem &p)> xform_fn;
	uint64_t cnt;

	if (this->type == t)
		return;

	xform_fn = this->getTransformFn(this->type, t);

	p0 = points_pop_begin();
	pe = points_pop_end();
	for(p=p0, cnt=0; p != pe; ++p){
		SpmCooElem p_new = SpmCooElem(*p);
		xform_fn(p_new);
		elems.push_back(p_new);
	}


	sort(elems.begin(), elems.end(), elem_cmp_less);
	SetRows(elems.begin(), elems.end());
	elems.clear();
	this->type = t;
}

std::ostream &operator<<(std::ostream &out, const SpmCooElem e)
{
	out << static_cast<CooElem>(e);
	if (e.pattern != NULL)
			out << "->[" << *(e.pattern) << "]";
	return out;
}

std::ostream &operator<<(std::ostream &out, const SpmRowElem &elem)
{
	out << elem.x;
	if (elem.pattern){
		out << *(elem.pattern);
	}
	return out;
}

std::ostream &operator<<(std::ostream &out, const SpmRowElems &elems)
{
	out << "row( ";
	BOOST_FOREACH(const SpmRowElem &elem, elems){
		out << elem << " [@" << &elem << "]  ";
	}
	out << ") @" << &elems ;
	return out;
}

std::ostream &operator<<(std::ostream &out, const SpmRows &rows)
{
	BOOST_FOREACH(const SpmRowElems &row, rows){
		out << row << std::endl;
	}
	return out;
}


#if 0
void SpmIdx::PrintRows(std::ostream &out)
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



void SpmIdx::genStats()
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

void SpmIdx::Print(std::ostream &out)
{
	SpmIdx::PointIter p, p_start, p_end;
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
	SpmIdx obj;
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
