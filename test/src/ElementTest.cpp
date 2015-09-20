/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/*
 * \file ElementTest.cpp
 * \brief Simple program for testing generic elements
 *
 * \author Vasileios Karakasis
 * \date 2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/CsxUtil.hpp>
#include <sparsex/internals/DynamicArray.hpp>
#include <sparsex/internals/Element.hpp>
#include <sparsex/internals/SparsePartition.hpp>
#include <sparsex/internals/Statistics.hpp>
#include <sparsex/internals/Utility.hpp>
#include <sparsex/internals/Xform.hpp>
#include <boost/container/vector.hpp>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;
using namespace sparsex::csx;
using namespace sparsex::utilities;

const size_t RleSize = 4;

int main()
{
    int row = 3;
    int col = 5;
    int delta = 2;
    vector<float> vals;
    for (size_t i = 0; i < RleSize; ++i)
        vals.push_back(5 + 2*i);

    //float *bogus = 0;
    Element<int, float> elem(row, col, vals.begin(), RleSize,
                              make_pair(Encoding::Horizontal, delta));
    Element<int, float> new_elem = elem;
    cout << elem << "\n";
    cout << new_elem << " (copied)\n";
    cout << GetPatternId(elem) << " (patt id)\n";
    assert(elem.IsPattern());

    // Test append in place
    boost::container::vector<Element<int, float> > elems;
    elems.emplace_back(row, col, vals.begin(), RleSize,
                       make_pair(Encoding::Horizontal, delta));
    cout << elems.back() << " (emplaced) \n";

    typename TransformFn<int>::type xform_fn =
        GetXformFn<int>(Encoding::Horizontal, Encoding::Diagonal);

    typename TransformFn<int>::type rev_xform_fn =
        GetXformFn<int>(Encoding::Diagonal, Encoding::Horizontal);

    Element<int, float> xformed = elem;
    xformed.Transform(xform_fn(elem.GetCoordinates(), 10, 10));
    cout << xformed << " (xformed)\n";

    Element<int, float> xformed_back = xformed;
    xformed_back.Transform(rev_xform_fn(xformed.GetCoordinates(),
                                                   10, 10));
    assert(xformed_back == elem);

    Element<int, float> single(1, 1, 4.3);
    cout << single << " (single)\n";
    xformed.GetMarker().Mark(PatternMarker::InPattern |
                             PatternMarker::PatternStart);
    cout << xformed.GetMarker().GetMark() << " (marked) \n";
    cout << xformed.GetMarker().IsMarked(PatternMarker::InPattern |
                                         PatternMarker::PatternStart)
         << " (marked?) \n";

    // Test insertion in dynamic array
    DynamicArray<Element<int, float> > da;
    da.Append(elem);
    da.Append(new_elem);
    da.Append(xformed);
    da.Append(xformed_back);
    da.Append(single);
    Element<int, float> *da_elems = da.TakeElems();
    da.GetAllocator().destroy(da_elems, da.GetSize());
    da.GetAllocator().deallocate(da_elems, da.GetSize());

    // Test statistics
    StatsCollection<StatsData> stats;
    stats.AppendStats(make_pair(Encoding::Horizontal, 1), StatsData(10, 2, 0));
    stats.AppendStats(make_pair(Encoding::Horizontal, 1), StatsData(10, 3, 0));
    stats.AppendStats(make_pair(Encoding::Horizontal, 2), StatsData(5, 3, 0));
    stats.AppendStats(make_pair(Encoding::Vertical, 1), StatsData(50, 2, 0));
    stats.AppendStats(make_pair(Encoding::Vertical, 0), StatsData(0, 0, 10));
//    stats.AppendStats(make_pair(Encoding::BlockRow2, 128), StatsData(768, 3, 0));
    // Test block splitter
    stats.AppendStats(make_pair(Encoding::BlockRow3, 2), StatsData(84, 14, 0));
    stats.AppendStats(make_pair(Encoding::BlockRow3, 3), StatsData(117, 13, 0));
    stats.AppendStats(make_pair(Encoding::BlockRow3, 4), StatsData(36, 3, 0));
    stats.AppendStats(make_pair(Encoding::BlockRow3, 5), StatsData(45, 3, 0));
    stats.AppendStats(make_pair(Encoding::BlockRow3, 6), StatsData(36, 2, 0));
    stats.AppendStats(make_pair(Encoding::BlockRow3, 12), StatsData(72, 2, 0));

    size_t nr_nzeros = 1000;
    set<Encoding::Instantiation> unfiltered;
    StatsDataScaler scaler(2.0);
    CoverageFilter cfilter(nr_nzeros, 0.1, unfiltered);
    BlockSplitter splitter(255, nr_nzeros, 0.1);
    // stats.ManipulateStats(scaler);
    stats.ManipulateStats(splitter);
    stats.ManipulateStats(cfilter);
    cout << stats.GetInstStats(make_pair(Encoding::Horizontal, 1)) << "\n";
    cout << stats.GetTypeStats(Encoding::Horizontal) << " (horiz)\n";
    cout << stats.GetTypeStats(Encoding::Vertical) << " (vert)\n";
    cout << stats.GetTypeStats(Encoding::BlockRow2) << " (br2)\n";
    cout << stats.GetTypeStats(Encoding::BlockRow3) << " (br3)\n";
    cout << stats.GetGlobalStats() << " (total)\n";
    cout << stats << "\n";
    cout << unfiltered.size() << " (total unfiltered inst)\n";
    //delete[] vals;

    // Test gcd and lcm
    cout << "gcd: " << gcd(44, 4) << "\n";
    cout << "abs: " << abs(-45.3) << "\n";
    cout << "lcm: " << lcm(44, 3) << "\n";
    cout << "iceil: " << iceil(5, 3) << "\n";

    // StatsCollection iterators
    StatsCollection<StatsData>::const_iterator si = stats.cbegin();
    StatsCollection<StatsData>::const_iterator se = stats.cend();
    cout << "Testing StatsCollection iterator\n";
    for (; si != se; ++si) {
        Encoding e((*si).first);
        cout << e.GetShortName() << ": " << (*si).second << "\n";
    }

    // Test append type stats
    StatsCollection<StatsData> new_stats;
    new_stats.AppendStats(Encoding::Horizontal, stats);
    cout << "new_stats: " << new_stats << "\n";

    stats.Clear();
    cout << stats.GetGlobalStats() << " (cleared)\n";
    return 0;
}
