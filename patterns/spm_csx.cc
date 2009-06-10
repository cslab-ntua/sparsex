#include <iterator>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

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

typedef struct {
	point_t point;
	drle_t *pattern; // if pattern is set, point is the starting point
} spm_elem_t;

static inline bool elem_cmp_larger(const spm_elem_t &e0, const spm_elem_t &e1)
{
	return (point_cmp_fn(e0.point, e1.point) > 1);
}

typedef std::vector<spm_elem_t> SpmIdxElems;
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

static inline void pnt_map_D(point_t &src, point_t &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	assert(nrows + src_x - src_y > 0);
	dst.y = nrows + src_x - src_y;
	dst.x = (src_x < src_y) ? src_x : src_y;
}
static inline void pnt_rmap_D(point_t &src, point_t &dst, uint64_t nrows)
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

static inline void pnt_map_rD(point_t &src, point_t &dst, uint64_t nrows)
{
	uint64_t src_x = src.x;
	uint64_t src_y = src.y;
	uint64_t dst_y;
	dst.y = dst_y = src_x + src_y - 1;
	dst.x = (dst_y <= nrows) ? src_x : src_x + nrows - dst_y;
}

static inline void pnt_rmap_rD(point_t &src, point_t &dst, uint64_t nrows)
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
	SpmIdxElems elems;

	SpmIdx() {type = NONE;};
	~SpmIdx() {};

	// load matrix from an MMF file
	void loadMMF(std::string mmf_file);
	void loadMMF(std::istream &in=std::cin);

	// Print Functions
	void Print(std::ostream &out=std::cout);
	void PrintRows(std::ostream &out=std::cout);

	void Transform(SpmIdxType type);
	void DeltaRLEncode(int limit=4, SpmIdxElems dstelems=elems);
	void Sort() {sort(elems.begin(), elems.end(), elem_cmp_larger);};
	void genStats();

	// Elements vector transformations:
	void XFormV() {
		FOREACH(spm_elem_t &e, this->elems){
			pnt_map_V(e.point, e.point);
		}
	};
	void XFormD() {
		FOREACH(spm_elem_t &e, this->elems){
			pnt_map_D(e.point, e.point, nrows);
		}
	};
	void XFormRD() {
		FOREACH(spm_elem_t &e, this->elems){
			pnt_map_rD(e.point, e.point, nrows);
		}
	};

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

void SpmIdx::loadMMF(std::istream &in)
{
	#define BUFFSIZE 1024
	char buff[BUFFSIZE];
	spm_elem_t elem;
	double val;
	int ret;
	SpmIdxElems elems;

	// header
	do {
		in.getline(buff, BUFFSIZE);
	} while (buff[0] == '#');
	ret = sscanf(buff, "%lu %lu %lu", &this->nrows, &this->ncols, &this->nnz);
	if (ret != 3){
		std::cerr << "mmf header error: sscanf" << std::endl;
		exit(1);
	}

	// body
	type = HORIZ;
	this->elems.resize(0);
	this->elems.reserve(this->nnz);
	elem.pattern = NULL;
	while (!in.getline(buff, BUFFSIZE).eof()){
		ret = sscanf(buff, "%lu %lu %lf", &elem.point.y, &elem.point.x, &val);
		if (ret != 3){
			std::cerr << "mmf body error: sscanf" << std::endl;
			exit(1);
		}
		this->elems.push_back(elem);
	}
	if (this->elems.size() != this->nnz){
		std::cerr << "mmf body error: size() != nnz" << std::endl;
		exit(1);
	}
	#undef BUFFSIZE
}

void SpmIdx::PrintRows(std::ostream &out)
{
	uint64_t prev_y=1;
	out << prev_y << ": ";
	FOREACH(spm_elem_t &e, this->elems){
		uint64_t y = e.point.y;
		uint64_t x = e.point.x;
		if (y != prev_y){
			out << std::endl << y << ": ";
			prev_y = y;
		}
		out << x << " ";
	}
	out << std::endl;
}

void SpmIdx::Print(std::ostream &out)
{
	FOREACH(spm_elem_t &e, this->elems){
		out << e.point << std::endl;
	}
}

void SpmIdx::Transform(SpmIdxType t)
{
	if (this->type == t)
		return;

	if (this->type != HORIZ)
		assert(false);

	switch(t) {
		case VERTICAL:
		XFormV();
		break;

		case DIAGONAL:
		XFormD();
		break;

		case REV_DIAGONAL:
		XFormRD();
		break;

		default:
		assert(false);
	}
	this->type = t;
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

int main(int argc, char **argv)
{
	//SpmIdx obj;
	//obj.loadMMF();
	//obj.Transform(DIAGONAL);
	//obj.Print();
}
