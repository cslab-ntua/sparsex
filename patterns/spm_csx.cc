#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <boost/function.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>

namespace bll = boost::lambda;

#define FOREACH BOOST_FOREACH

namespace csx {

typedef struct {
	uint64_t x, y;
} point_t;

static inline int point_cmp_fn(const point_t &p0, const point_t &p1)
{
	int64_t ret;
	ret = p0.y - p1.y;
	if (ret == 0){
		ret = p0.x - p1.x;
	}
	if (ret > 0){
		return 1;
	} else if (ret < 0){
		return -1;
	} else {
		return 0;
	}
}

typedef struct {
	uint32_t size, drle_len;
	enum {DRLE_H=1, DRLE_V, DRLE_D, DRLE_RD} type;
} drle_t;


typedef struct : public point_t {
	drle_t *pattern;
} spm_coo_elem_t;

typedef struct {
	uint64_t x;
	drle_t *pattern;
} spm_row_elem_t;

static inline bool elem_cmp_less(const spm_coo_elem_t &e0,
                                   const spm_coo_elem_t &e1)
{
	int ret;
	ret = point_cmp_fn(static_cast<point_t>(e0), static_cast<point_t>(e1));
	return (ret < 0);
}

typedef std::vector<point_t> SpmPoints;
typedef std::iterator<std::forward_iterator_tag, point_t> SpmPointIter;
typedef std::vector<spm_coo_elem_t> SpmCooElems;
typedef std::vector<spm_row_elem_t> SpmRowElems;
typedef std::vector<SpmRowElems> SpmRows;

typedef enum {NONE=0, HORIZ, VERTICAL, DIAGONAL, REV_DIAGONAL} SpmIdxType;
typedef boost::function<void (point_t &p)> TransformFn;

// mappings for vertical transformation
static inline void pnt_map_V(point_t &src, point_t &dst)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	dst.x = src_y;
	dst.y = src_x;
}
static inline void pnt_rmap_V(point_t &src, point_t &dst)
{
	pnt_map_V(src, dst);
}

static inline void pnt_map_D(const point_t &src, point_t &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	assert(nrows + src_x - src_y > 0);
	dst.y = nrows + src_x - src_y;
	dst.x = (src_x < src_y) ? src_x : src_y;
}
static inline void pnt_rmap_D(const point_t &src, point_t &dst, uint64_t nrows)
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

static inline void pnt_map_rD(const point_t &src, point_t &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	uint64_t dst_y;
	dst.y = dst_y = src_x + src_y - 1;
	dst.x = (dst_y <= nrows) ? src_x : src_x + nrows - dst_y;
}

static inline void pnt_rmap_rD(const point_t &src, point_t &dst, uint64_t nrows)
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


class SpmIdx {
public:
	uint64_t nrows, ncols, nnz;
	SpmIdxType type;
	SpmRows rows;

	SpmIdx() {type = NONE;};
	~SpmIdx() {};

	// load matrix from an MMF file
	void loadMMF(std::string mmf_file);
	void loadMMF(std::istream &in=std::cin);

	// Print Functions
	void Print(std::ostream &out=std::cout);
	void PrintRows(std::ostream &out=std::cout);

	template <typename IterT>
	void SetRows(IterT pnts_start, IterT pnts_end);

	inline TransformFn getTransformFn(SpmIdxType type);
	inline TransformFn getRTransformFn(SpmIdxType type);
	void Transform(SpmIdxType type);

	// iterators that return a point_t item
	class PointIter;
	PointIter points_begin();
	PointIter points_end();

	// same with PointIter, but removes elements
	class PointPoper;
	PointPoper points_pop_begin();
	PointPoper points_pop_end();

	#if 0
	void DeltaRLEncode(int limit=4, SpmIdxElems dstelems=elems);
	void Sort() {sort(elems.begin(), elems.end(), elem_cmp_larger);};
	void genStats();

	void XFormV() {
		FOREACH(spm_coo_elem_t &e, this->elems){
			pnt_map_V(e.point, e.point);
		}
	};
	void XFormD() {
		FOREACH(spm_coo_elem_t &e, this->elems){
			pnt_map_D(e.point, e.point, nrows);
		}
	};
	void XFormRD() {
		FOREACH(spm_coo_elem_t &e, this->elems){
			pnt_map_rD(e.point, e.point, nrows);
		}
	};
	#endif

};


class SpmIdx::PointIter
: public std::iterator<std::forward_iterator_tag, point_t>
{
public:
	uint64_t row_idx;
	uint64_t elm_idx;
	SpmIdx   *spm;
	PointIter(): row_idx(0), elm_idx(0), spm(NULL) {;}
	PointIter(uint64_t ridx, uint64_t eidx, SpmIdx *s)
	: row_idx(ridx), elm_idx(eidx), spm(s) {
		// find the first non zero row
		while (row_idx < spm->nrows && spm->rows[row_idx].size() == 0) {
			row_idx++;
		}
	}

	bool operator==(const PointIter &x){
		return (spm == x.spm)
			   && (row_idx == x.row_idx)
			   && (elm_idx == x.elm_idx);
	}

	bool operator!=(const PointIter &x){
		return !(*this == x);
	}

	void operator++() {
		uint64_t rows_nr, row_elems_nr;

		rows_nr = spm->nrows;
		row_elems_nr = spm->rows[row_idx].size();

		// equality means somebody did a ++ on an ended iterator
		assert(row_idx < rows_nr);

		elm_idx++;
		assert(elm_idx <= row_elems_nr);
		if (elm_idx < row_elems_nr){
			return ;
		}

		// change row
		do {
			elm_idx = 0;
			row_idx++;
		} while (row_idx < rows_nr && spm->rows[row_idx].size() == 0);
	}

	spm_coo_elem_t operator*(){
		spm_coo_elem_t ret;
		ret.y = row_idx + 1;
		ret.x = spm->rows[row_idx][elm_idx].x;
		ret.pattern = spm->rows[row_idx][elm_idx].pattern;
		return ret;
	}
};

class SpmIdx::PointPoper : public SpmIdx::PointIter
{
public:
	PointPoper() : PointIter() { ; }
	PointPoper(uint64_t ridx, uint64_t eidx, SpmIdx *s)
	: PointIter(ridx, eidx, s) { ; }

	void operator++() {
		uint64_t rows_nr, row_elems_nr;

		rows_nr = spm->nrows;
		row_elems_nr = spm->rows[row_idx].size();

		// equality means somebody did a ++ on an ended iterator
		assert(row_idx < rows_nr);

		elm_idx++;
		assert(elm_idx <= row_elems_nr);
		if (elm_idx < row_elems_nr){
			return ;
		}

		// remove all elements from previous row
		spm->rows[row_idx].resize(0);
		// change row
		do {
			elm_idx = 0;
			row_idx++;
		} while (row_idx < rows_nr && spm->rows[row_idx].size() == 0);
	}
};


} // end of csx namespace

using namespace csx;

std::istream &operator>>(std::istream &in, SpmIdx &obj)
{
	obj.loadMMF(in);
	return in;
}

std::ostream &operator<<(std::ostream &out, point_t p)
{
	out << "(" << std::setw(2) << p.y << "," << std::setw(2) << p.x << ")";
	return out;
}

std::ostream &operator<<(std::ostream &out, SpmIdx::PointIter pi)
{
	out << "<" << std::setw(2) << pi.row_idx << "," << std::setw(2) << pi.elm_idx << ">";
	return out;
}

SpmIdx::PointIter SpmIdx::points_begin()
{
	return PointIter(0, 0, this);
}

SpmIdx::PointIter SpmIdx::points_end()
{
	return PointIter(nrows, 0, this);
}

SpmIdx::PointPoper SpmIdx::points_pop_begin()
{
	return PointPoper(0, 0, this);
}

SpmIdx::PointPoper SpmIdx::points_pop_end()
{
	return PointPoper(nrows, 0, this);
}

static inline spm_row_elem_t mk_row_elem(const point_t &p)
{
	static spm_row_elem_t ret;
	ret.x = p.x;
	ret.pattern = NULL;
	return ret;
}

static inline spm_row_elem_t mk_row_elem(const spm_coo_elem_t &p)
{
	static spm_row_elem_t ret;
	ret.x = p.x;
	ret.pattern = p.pattern;
	return ret;
}

template <typename IterT>
void SpmIdx::SetRows(IterT pnts_start, IterT pnts_end)
{
	SpmRowElems *row_elems;
	spm_row_elem_t elem;
	uint64_t row_prev;

	//this->rows.resize(this->nrows);
	this->rows.resize(1);
	row_elems = &this->rows.at(0);
	row_prev = 1;

	// XXX
	IterT point;
	for (point = pnts_start; point < pnts_end; ++point){
		uint64_t row = point->y;
		if (row != row_prev){
			assert(row > this->rows.size());
			this->rows.resize(row);
			row_elems = &this->rows.at(row - 1);
			row_prev = row;
		}
		elem = mk_row_elem(*point);
		row_elems->push_back(elem);
	}
}

void SpmIdx::loadMMF(std::istream &in)
{
	char buff[1024];
	double val;
	uint64_t cnt;
	int ret;
	point_t *point;
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
	type = HORIZ;
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

inline TransformFn SpmIdx::getRTransformFn(SpmIdxType type)
{
	boost::function<void (point_t &p)> ret;
	ret = NULL;
	switch(type) {
		case HORIZ:
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

inline TransformFn SpmIdx::getTransformFn(SpmIdxType type)
{
	boost::function<void (point_t &p)> ret;
	ret = NULL;
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

		default:
		assert(false);
	}
	return ret;
}

void SpmIdx::Transform(SpmIdxType t)
{
	PointPoper p0, pe, p;
	SpmCooElems elems;
	boost::function<void (point_t &p)> xform_fn, rxform_fn;
	uint64_t cnt;

	if (this->type == t)
		return;

	xform_fn = getTransformFn(t);
	rxform_fn = getRTransformFn(this->type);
	if (rxform_fn != NULL){
		xform_fn = bll::bind(xform_fn, (bll::bind(rxform_fn, bll::_1), bll::_1));
	}

	p0 = points_pop_begin();
	pe = points_pop_end();
	elems.resize(this->nnz);
	for(p=p0, cnt=0; (p != pe) && (cnt < this->nnz); ++p){
		spm_coo_elem_t *myp = &elems[cnt++];
		*myp = *p;
		xform_fn(*myp);
	}
	assert(cnt == this->nnz);

	sort(elems.begin(), elems.end(), elem_cmp_less);
	SetRows(elems.begin(), elems.end());
	elems.clear();
	this->type = t;
}

#if 0
void SpmIdx::PrintRows(std::ostream &out)
{
	uint64_t prev_y=1;
	out << prev_y << ": ";
	FOREACH(spm_coo_elem_t &e, this->elems){
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

void SpmIdx::DeltaRLEncode(int limit, SpmIdxElems dstelems)
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

int main(int argc, char **argv)
{
	SpmIdx obj;
	obj.loadMMF();
	obj.Transform(DIAGONAL);

	//obj.Print();
	//obj.Print();
}
