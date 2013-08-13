/* -*- C++ -*-
 *
 * rcm.h -- Reverse Cuthill-Mckee Ordering Interface
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_RCM_H__
#define CSX_RCM_H__

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "mmf.h"
#include "csr.h"
#include "csr_iterator.h"
#include "SparseInternal.h"

using namespace boost;

typedef pair<size_t, size_t> Pair;
typedef adjacency_list<vecS, vecS, undirectedS, 
                       property<vertex_color_t, default_color_type,
                       property<vertex_degree_t,int> > > Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertices_size_type size_type;

/**
 *  @return perm      permutation from the old ordering to the new one.
 *  @return inv_perm  permutation from the new ordering to the old one.
 *  @param  graph     the incidence graph of the matrix.
 */
void FindPerm(vector<size_t>& perm, vector<size_t>& inv_perm, Graph& graph);

/**
 *  Loads matrix from an MMF file and applies the reverse Cuthill-Mckee 
 *  reordering algorithm
 *
 *  @param file_name   name of the file where the matrix is kept.
 *  @param rt_config   runtime configuration.
 *  @param mat         handler of MMF class.
 *  @param graph       the incidence graph of the matrix.
 *  @param perm        permutation from the old ordering to the new one.
 *  @return            spm class object with the characteristics of the matrix.
 */
template<typename IndexType, typename ValueType>
void ReorderMat_MMF(MMF<IndexType, ValueType>& mat, const vector<size_t>& perm);
template<typename IndexType, typename ValueType>
Graph& ConstructGraph_MMF(Graph& graph, MMF<IndexType, ValueType>& mat);
template<typename IndexType, typename ValueType>
void DoReorder_RCM(MMF<IndexType, ValueType>& mat, vector<size_t>& inv_perm);

template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadMMF_RCM(const char *filename);
template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadMMF_RCM_mt(const char *file_name,
                                                     const long nr_threads = 0);

/**
 *  Loads matrix from CSR format and applies the reverse Cuthill-Mckee 
 *  reordering algorithm
 *
 *  @param rowptr      array "rowptr" of CSR format.
 *  @param colind      array "colind" of CSR format.
 *  @param values      array "values" of CSR format.
 *  @param nr_rows     number of rows.
 *  @param nr_cols     number of columns.
 *  @param zero_based  indexing.
 *  @param rt_config   runtime configuration.
 *  @param perm        permutation from the old ordering to the new one.
 *  @param inv_perm    permutation from the new ordering to the old one.
 *  @param graph       the incidence graph of the matrix.
 *  @return            spm class object with the characteristics of the matrix.
 */
template<typename IndexType, typename ValueType>
void ReorderMat_CSR(CSR<IndexType, ValueType>& mat, const vector<size_t>& perm,
                    vector<size_t>& inv_perm);
template<typename IterT>
Graph& ConstructGraph_CSR(Graph& graph, IterT& iter, const IterT& iter_end,
                          size_t nr_nzeros, bool symmetric);
template<typename IndexType, typename ValueType>
void DoReorder_RCM(CSR<IndexType, ValueType>& mat, vector<size_t>& inv_perm);

template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadCSR_RCM_mt(IndexType *rowptr,
                                                     IndexType *colind,
                                                     ValueType *values,
                                                     IndexType nr_rows,
                                                     IndexType nr_cols,
                                                     bool zero_based,
                                                     const long nr_threads = 0);


/*
 * Implementation of RCM interface
 */
void FindPerm(vector<size_t>& perm, vector<size_t>& invperm, Graph& graph)
{
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
    invperm.reserve(num_vertices(graph));

    // Reverse Cuthill Mckee Ordering
    cuthill_mckee_ordering(graph, inv_perm.rbegin(), get(vertex_color, graph),
                           make_degree_map(graph));

    for (size_t i = 0; i != inv_perm.size(); ++i) {
        perm[index_map[inv_perm[i]]] = i;
        invperm.push_back(inv_perm[i]);
    }

    // cout << "Permutation: ";
    // for (size_t i = 0; i < inv_perm.size(); i++) {
    //     cout << inv_perm[i] << " ";
    // }
    // cout << endl;

    size_t nb = bandwidth(graph, make_iterator_property_map
                          (&perm[0], index_map, perm[0]));
                          cout << "Final Bandwidth: " << nb << endl;
}

template<typename IndexType, typename ValueType>
Graph& ConstructGraph_MMF(Graph& graph, MMF<IndexType, ValueType>& mat)
{
    // make a pessimistic guess for nr_edges
    size_t nr_edges;
    if (mat.IsSymmetric()) {    //if main diagonal is full (nr_nzeros - nr_rows) / 2
        nr_edges = mat.GetNrNonzeros() / 2;
    } else {
        nr_edges = mat.GetNrNonzeros();
    }
        
    Pair *edges = (Pair *) malloc(nr_edges * sizeof(Pair));

    typename MMF<IndexType, ValueType>::iterator iter = mat.begin();
    typename MMF<IndexType, ValueType>::iterator iter_end = mat.end();
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

template<typename IndexType, typename ValueType>
void ReorderMat_MMF(MMF<IndexType, ValueType>& mat, const vector<size_t>& perm)
{    
    IndexType row, col;

    for (IndexType i = 0; i < mat.GetNrNonzeros(); i++) {
        mat.GetCoordinates(i, row, col);
        mat.SetCoordinates(i, perm[row-1] + 1, perm[col-1] + 1);
    }
    mat.Sort();
    //mat.Print();
}

template<typename IndexType, typename ValueType>
void DoReorder_RCM(MMF<IndexType, ValueType>& mat, vector<size_t> &perm)
{    
    vector<size_t> inv_perm;

    cout << "=== START REORDERING ===\n";
    // Construct graph
    Graph graph(mat.GetNrRows());
    graph = ConstructGraph_MMF(graph, mat);
    // Find permutation
    FindPerm(perm, inv_perm, graph);  
    // Reorder original matrix
    ReorderMat_MMF(mat, perm);

    cout << "=== END REORDERING ===\n";
}

template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadMMF_RCM_mt(const char *mmf_file,
                                                     const long nr_threads)
{
    SparseInternal<IndexType, ValueType> *ret = 0;
    MMF<IndexType, ValueType> mat(mmf_file);
    vector<size_t> inv_perm;

    // cout << "=== START REORDERING ===\n";
    // Reorder matrix
    DoReorder_RCM(mat, inv_perm);
    // Create internal representation
    ret = SparseInternal<IndexType, ValueType>::DoLoadMatrix(mat, nr_threads);
    // ret->Print(std::cout);
    // cout << "=== END REORDERING ===\n";
    return ret;
}

template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadMMF_RCM(const char *mmf_file)
{
    SparseInternal<IndexType, ValueType> *ret;
    ret = LoadMMF_RCM_mt<IndexType, ValueType>(mmf_file, 1);
    return ret;
}

template<typename IterT>
Graph& ConstructGraph_CSR(Graph& graph, IterT& iter, const IterT& iter_end,
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
    if (symmetric) {
        for (;iter != iter_end; ++iter) {
            if ((*iter).row < (*iter).col) {
                edges[index++] = Pair(((*iter).row) - 1, ((*iter).col) - 1);
            }
        }
    } else {
        for (;iter != iter_end; ++iter) {
            if ((*iter).row != (*iter).col) {
                edges[index++] = Pair(((*iter).row) - 1, ((*iter).col) - 1);
            }
        }
    }
    // index -> actual nr_edges
    assert(index <= nr_edges);

    if (index == 0) {
        cerr << "No reordering available for this matrix: "
             << "all non-zero elements on main diagonal" << endl;
        exit(1);
    }

    for (size_t i = 0; i < index; i++) {
        add_edge(edges[i].first, edges[i].second, graph);
    }
    free(edges);

    return graph;
}

template<typename IndexType, typename ValueType>
void ReorderMat_CSR(CSR<IndexType, ValueType>& mat, const vector<size_t>& perm,
                    vector<size_t>& inv_perm)
{
    // Apply permutation only to colind
    for (IndexType i = 0; i < mat.GetNrNonzeros(); i++) {
        mat.colind_[i] = perm[mat.colind_[i] - !mat.IsZeroBased()] +
            !mat.IsZeroBased();
    }

    // Sort colind and values per row
    IndexType row_start, length;
    for (IndexType i = 0; i < mat.GetNrRows(); i++) {
        row_start = mat.rowptr_[i] - !mat.IsZeroBased();
        length = mat.rowptr_[i+1] - mat.rowptr_[i];
        std::sort(get_CSR_iterator<IndexType*, ValueType*>
                  (&mat.colind_[row_start], &mat.values_[row_start]), 
                  get_CSR_iterator<IndexType*, ValueType*>
                  (&mat.colind_[row_start] + length,
                   &mat.values_[row_start] + length), 
                  CSR_Comp<IndexType*, ValueType*>());
    }
    mat.SetReordered(inv_perm);
    //assert(inv_perm.capacity() == 0);
}

template<typename IndexType, typename ValueType>
void DoReorder_RCM(CSR<IndexType, ValueType>& mat, vector<size_t>& perm)
{
    typename CSR<IndexType, ValueType>::iterator iter = mat.begin();
    typename CSR<IndexType, ValueType>::iterator iter_end = mat.end();
    vector<size_t> inv_perm;

    cout << "=== START REORDERING ===\n";
    // Construct graph
    Graph graph(mat.GetNrRows());
    graph = ConstructGraph_CSR(graph, iter, iter_end, mat.GetNrNonzeros(),
                               mat.IsSymmetric());
    // Find permutation
    FindPerm(perm, inv_perm, graph);
    // Reorder original matrix
    ReorderMat_CSR<IndexType, ValueType>(mat, perm, inv_perm);
    cout << "=== END REORDERING ===\n";
}

template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadCSR_RCM_mt(IndexType *rowptr,
                                                     IndexType *colind,
                                                     ValueType *values,
                                                     IndexType nr_rows,
                                                     IndexType nr_cols,
                                                     bool zero_based,
                                                     const long nr_threads)
{
    SparseInternal<IndexType, ValueType> *spi = NULL;
    CSR<IndexType, ValueType> mat(rowptr, colind, values, 
                                  nr_rows, nr_cols, zero_based);
    vector<size_t> inv_perm;

    // cout << "=== START REORDERING ===\n";

    // Reorder original matrix
    DoReorder_RCM(mat, inv_perm);
    // Create internal representation
    spi = SparseInternal<IndexType, ValueType>::DoLoadMatrix(mat, nr_threads);
    //spms->Print();

    // cout << "=== END REORDERING ===\n";

    return spi;
}

#endif // RCM_H__
