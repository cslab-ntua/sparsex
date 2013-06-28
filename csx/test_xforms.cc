/*
 * test_xforms.cc --  Testing transformations.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "sparsepartition.h"

#define FOREACH BOOST_FOREACH

namespace bll = boost::lambda;
using namespace csx;

template<typename IndexType, typename ValueType>
void PrintTransform(IndexType y, IndexType x,
                    typename TransformFnType<IndexType, ValueType>::
                    TransformFn xform_fn,
                    std::ostream &out)
{
    IndexType i,j;
    CooElem<IndexType, ValueType> p;

    for (i = 1; i <= y; i++) {
        for (j = 1; j <= x;  j++) {
            p.row = i;
            p.col = j;
            xform_fn(p);
            out << p << " ";
        }

        out << std::endl;
    }
}

template<typename IndexType, typename ValueType>
void PrintDiagTransform(IndexType y, IndexType x, std::ostream &out)
{
    boost::function<void (CooElem<IndexType, ValueType> &p)> xform_fn;

    xform_fn = bll::bind(pnt_map_D<IndexType, ValueType>, bll::_1, bll::_1, y);
    PrintTransform(y, x, xform_fn, out);
}

template<typename IndexType, typename ValueType>
void PrintRDiagTransform(IndexType y, IndexType x, std::ostream &out)
{
    boost::function<void (CooElem<IndexType, ValueType> &p)> xform_fn;

    xform_fn = bll::bind(pnt_map_rD<IndexType, ValueType>, bll::_1, bll::_1, x);
    PrintTransform(y, x, xform_fn, out);
}

template<typename IndexType, typename ValueType>
void TestTransform(IndexType y, IndexType x,
                   typename TransformFnType<IndexType, ValueType>::
                   TransformFn xform_fn,
                   typename TransformFnType<IndexType, ValueType>::
                   TransformFn rxform_fn)
{
    IndexType i,j;
    CooElem<IndexType, ValueType> p0, p1;
    
    for (i = 1; i <= y; i++) {
        for (j = 1; j <= x;  j++){
            p0.row = i;
            p0.col = j;
            xform_fn(p0);
            p1.row = p0.row;
            p1.col = p0.col;
            rxform_fn(p1);
            if (p1.row != i || p1.col != j) {
                cout << "Error for " << i << "," << j << endl;
                cout << "Transformed: " << p0.row << "," << p0.col << endl;
                cout << "rTransformed:" << p1.row << "," << p1.col << endl;
                exit(1);
            }
        }
    }
}

template<typename IndexType, typename ValueType>
void TestRDiagTransform(IndexType y, IndexType x)
{
    boost::function<void (CooElem<IndexType, ValueType> &p)> xform_fn, rxform_fn;
    
    xform_fn = bll::bind(pnt_map_rD<IndexType, ValueType>, bll::_1, bll::_1, x);
    rxform_fn = bll::bind(pnt_rmap_rD<IndexType, ValueType>,
                          bll::_1, bll::_1, x);
    TestTransform(y, x, xform_fn, rxform_fn);
}

#if 0

int main(int argc, char **argv)
{
	PrintRDiagTransform(5, 5, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(10, 5, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(5, 10, std::cout);
	std::cout << std::endl;

	TestRDiagTransform(5,5);
	TestRDiagTransform(10,5);
	TestRDiagTransform(5,10);
}

void printXforms()
{
	cout << "Diagonal" << endl;
	PrintDiagTransform(10, 10, cout);
	cout << endl;
	PrintDiagTransform(5, 10, cout);
	cout << endl;
	PrintDiagTransform(10, 5, cout);
	cout << endl;
	cout << "Reverse Diagonal" << endl;
	PrintRDiagTransform(10, 10, cout);
	cout << endl;
	PrintRDiagTransform(5, 10, cout);
	cout << endl;
	PrintRDiagTransform(10, 5, cout);
	cout << endl;
}

template<typename IndexType, typename ValueType>
void TestDiagTransform(IndexType y, IndexType x)
{
	boost::function<void (point_t &p)> xform_fn, rxform_fn;
	xform_fn = boost::bind(pnt_map_D<IndexType, ValueType>, _1, _1, y);
	rxform_fn = boost::bind(pnt_rmap_D<IndexType, ValueType>, _1, _1, y);
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

void test_drle()
{
	vector<int> v_in, delta;
	vector< RLE<int> > drle;

	v_in.resize(16);
	for (int i=0; i<16; i++){
		v_in[i] = i;
	}
	FOREACH(int &v, v_in){
		cout << v << " ";
	}
	cout << endl;

	delta = DeltaEncode(v_in);
	FOREACH(int &v, delta){
		cout << v << " ";
	}
	cout << endl;

	drle = RLEncode(delta);
	FOREACH(RLE<int> &v, drle){
		cout << "(" << v.val << "," << v.freq << ") ";
	}
	cout << endl;
}

int main(int argc, char **argv)
{
	SparsePartition<uint64_t, double> obj;
	obj.LoadMMF();
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

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
