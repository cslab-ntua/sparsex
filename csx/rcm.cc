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

#include "runtime.h"
#include "mmf.h"
#include "csr.h"
#include "csr_iterator.h"

using namespace std;
using namespace boost;

typedef pair<size_t, size_t> Pair;
typedef adjacency_list<vecS, vecS, undirectedS, 
                       property<vertex_color_t, default_color_type,
                       property<vertex_degree_t,int> > > Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertices_size_type size_type;

void FindPerm(vector<size_t> &perm, vector<size_t> &invperm, Graph &graph)
{
    cout << "=== RCM REORDERING ===\n";
    graph_traits<Graph>::vertex_iterator ui, ui_end;
    property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, graph);
    for (boost::tie(ui, ui_end) = vertices(graph); ui != ui_end; ++ui) {
        deg[*ui] = degree(*ui, graph);
    }

    // Original ordering
    property_map<Graph, vertex_index_t>::type 
        index_map = get(vertex_index, graph);

    size_t ob = bandwidth(graph);
    cout << "Original Bandwidth: " << ob << endl;

    vector<Vertex> inv_perm(num_vertices(graph));
    perm.resize(num_vertices(graph));
    invperm.resize(num_vertices(graph));

    // Reverse Cuthill Mckee Ordering
    cuthill_mckee_ordering(graph, inv_perm.rbegin(), get(vertex_color, graph),
                           make_degree_map(graph));

    for (size_t i = 0; i != inv_perm.size(); ++i) {
        perm[index_map[inv_perm[i]]] = i;
        invperm[i] = inv_perm[i];
    }

    cout << "Permutation: ";
    for (size_t i = 0; i < inv_perm.size(); i++) {
        cout << inv_perm[i] << " ";
    }
    cout << endl;

    size_t nb = bandwidth(graph, make_iterator_property_map
                          (&perm[0], index_map, perm[0]));
                          cout << "Final Bandwidth: " << nb << endl;
}

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
            if ((*iter).y < (*iter).x)
                edges[index++] = Pair(((*iter).y) - 1, ((*iter).x) - 1);
        }
    } else {
        for (;iter != iter_end; ++iter)
            if ((*iter).y < (*iter).x)
                edges[index++] = Pair(((*iter).y) - 1, ((*iter).x) - 1);
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

SPM *LoadMMF_RCM(const char *file_name, RuntimeContext &rt_config)
{
    SPM *spms = NULL;
    ifstream in(file_name);
    MMF mat(in);
    vector<size_t> perm, inv_perm;

    // Construct graph
    Graph G;
    G = ConstructGraphMMF(G, mat);
    // Find permutation
    FindPerm(perm, inv_perm, G);  
    // Reorder original matrix
    ReorderMatMMF(mat, perm);
    // Create internal representation
    spms = SPM::DoLoadMatrix<MMF>(mat, rt_config.GetNrThreads());
    //spms->Print();
    in.close();

    return spms;
}

template<typename IterT>
Graph& ConstructGraphCSR(Graph &G, IterT &iter, const IterT &iter_end,
                         size_t nr_nzeros, bool symmetric)
{
    // make a pessimistic guess for nr_edges
    size_t nr_edges;
    if (symmetric) {    //if main diagonal is full (nr_nzeros - nr_rows) / 2
        nr_edges = nr_nzeros / 2;
    } else {
        nr_edges = nr_nzeros;
    }
        
    Pair *edges = (Pair *) malloc(nr_edges * sizeof(Pair));

    size_t index = 0;
    for (;iter != iter_end; ++iter) {
        if ((*iter).y < (*iter).x) {
                edges[index++] = Pair(((*iter).y) - 1, ((*iter).x) - 1);
        }
    }
    // index -> actual nr_edges
    assert(index <= nr_edges);
    cout << index << "edges\n";
    if (index == 0) {
        cerr << "No reordering available for this matrix: "
             << "all non-zero elements on main diagonal" << endl;
        exit(1);
    }

    for (size_t i = 0; i < index; ++i) {
        add_edge(edges[i].first, edges[i].second, G);
    }
    free(edges);

    return G;
}

template<typename IndexType, typename ValueType>
void ReorderMatCSR(CSR<IndexType, ValueType> &mat, vector<size_t> &perm,
                   vector<size_t> &inv_perm)
{
    IndexType *rowptr;
    rowptr = (IndexType *) malloc(sizeof(IndexType) * (mat.GetNrRows()+1));
    IndexType *colind;
    colind = (IndexType *)malloc(sizeof(IndexType) * (mat.GetNrNonzeros()));
    ValueType *values;
    values = (ValueType *)malloc(sizeof(ValueType) * (mat.GetNrNonzeros()));

    IndexType index = 0;
    for (size_t i = 0; i < mat.GetNrRows(); i++) {
        rowptr[i] = index + !mat.IsZeroBased();
        size_t k = inv_perm[i];
        IndexType row_start = index;
        IndexType length = mat.rowptr_[k+1] - mat.rowptr_[k];
        for (size_t j = mat.rowptr_[k]; j < mat.rowptr_[k+1]; j++) {
            values[index] = mat.values_[j-!mat.IsZeroBased()];
            colind[index++] = perm[mat.colind_[j-!mat.IsZeroBased()] - 
                                       !mat.IsZeroBased()] + !mat.IsZeroBased();
            //cout << new_colind[index-1] << " " << new_values[index-1] << endl;
        }   
        // Sort row
        std::sort(get_CSR_iterator<IndexType*, ValueType*>
                  (&colind[row_start], &values[row_start]), 
                  get_CSR_iterator<IndexType*, ValueType*>
                  (&colind[row_start] + length, &values[row_start] + length), 
                  CSR_Comp<IndexType*, ValueType*>());
    }
    rowptr[mat.GetNrRows()] = index + !mat.IsZeroBased();
    // FIXME: memory leak
    mat.rowptr_ = rowptr;
    mat.colind_ = colind;
    mat.values_ = values;
    rowptr = NULL;
    colind = NULL;
    values = NULL;
}

template<typename IndexType, typename ValueType>
SPM *LoadCSR_RCM(IndexType *rowptr, IndexType *colind, ValueType *values,
                 size_t nr_rows, size_t nr_cols, bool zero_based,
                 RuntimeContext &rt_config)
{
    SPM *spms = NULL;
    CSR<IndexType, ValueType> mat(rowptr, colind, values, 
                                  nr_rows, nr_cols, zero_based);

    typename CSR<IndexType, ValueType>::iterator iter = mat.begin();
    typename CSR<IndexType, ValueType>::iterator iter_end = mat.end();
    vector<size_t> perm, inv_perm;

    // Construct graph
    Graph G(nr_rows);
    G = ConstructGraphCSR(G, iter, iter_end, mat.GetNrNonzeros(),
                          mat.IsSymmetric());
    // Find permutation
    FindPerm(perm, inv_perm, G);
    // Reorder original matrix
    ReorderMatCSR<IndexType, ValueType>(mat, perm, inv_perm);
    // Create internal representation
    spms = SPM::DoLoadMatrix<CSR<IndexType, ValueType> >
        (mat, rt_config.GetNrThreads());

    return spms;
}


