/* -*- C++ -*-
 * 
 * rcm.cc -- Reverse Cuthill-McKee Ordering Interface
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "rcm.h"

using namespace std;
using namespace boost;

Graph& ConstructGraphMMF(Graph &graph, MMF &mat)
{
    // make a pessimistic guess for nr_edges
    size_t nr_edges;
    if (mat.IsSymmetric()) {    //if main diagonal is full (nr_nzeros - nr_rows) / 2
        nr_edges = mat.GetNrNonzeros() / 2;
    } else {
        nr_edges = mat.GetNrNonzeros();
    }
        
    Pair *edges = (Pair *) malloc(nr_edges * sizeof(Pair));

    MMF::iterator iter = mat.begin();
    MMF::iterator iter_end = mat.end();
    size_t index = 0;

    if (!mat.IsSymmetric() && !mat.IsColWise()) {
        mat.InitMatrix(mat.GetNrNonzeros());
        for (;iter != iter_end; ++iter) {
            mat.InsertElement(*iter);
            if ((*iter).row != (*iter).col)
                edges[index++] = Pair(((*iter).row) - 1, ((*iter).col) - 1);
        }
    } else {
        for (;iter != iter_end; ++iter)
            if ((*iter).row < (*iter).col)
                edges[index++] = Pair(((*iter).row) - 1, ((*iter).col) - 1);
    }
    // index -> actual nr_edges
    assert(index <= nr_edges);

    if (index == 0) {
        cerr << "No reordering available for this matrix: "
             << "all non-zero elements on main diagonal" << endl;
        exit(1);
    }

    if (!mat.IsSymmetric() && !mat.IsColWise()) 
        mat.SetReordered();

    for (size_t i = 0; i < index; ++i) {
        add_edge(edges[i].first, edges[i].second, graph);
    }
    free(edges);

    return graph;
}

void ReorderMatMMF(MMF &mat, const vector<size_t> &perm)
{    
    uint64_t row, col;

    for (size_t i = 0; i < mat.GetNrNonzeros(); i++) {
        mat.GetCoordinates(i, row, col);
        mat.SetCoordinates(i, perm[row-1] + 1, perm[col-1] + 1);
    }
    mat.Sort();
    //mat.Print();
}

SPM *LoadMMF_RCM(const char *file_name, 
                 const RuntimeContext &rt_config)
{
    SPM *spms = NULL;
    ifstream in(file_name);
    MMF mat(in);
    vector<size_t> perm, inv_perm;

    cout << "=== START REORDERING ===\n";

    // Construct graph
    Graph graph(mat.GetNrRows());
    graph = ConstructGraphMMF(graph, mat);
    // Find permutation
    FindPerm(perm, inv_perm, graph);  
    // Reorder original matrix
    ReorderMatMMF(mat, perm);
    // Create internal representation
    spms = SPM::DoLoadMatrix<MMF>(mat, rt_config.GetNrThreads());
    //spms->Print();
    in.close();

    cout << "=== END REORDERING ===\n";

    return spms;
}


