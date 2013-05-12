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

#ifndef RCM_H__
#define RCM_H__

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "runtime.h"
#include "mmf.h"
#include "csr.h"
#include "rcm_impl.h"

using namespace boost;

typedef pair<size_t, size_t> Pair;
typedef adjacency_list<vecS, vecS, undirectedS, 
                       property<vertex_color_t, default_color_type,
                       property<vertex_degree_t,int> > > Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertices_size_type size_type;

// FIXME: return permutation in order to reorder the rest of the linear system

/**
 *  @return perm      permutation from the old ordering to the new one.
 *  @return inv_perm  permutation from the new ordering to the old one.
 *  @param  graph     the incidence graph of the matrix.
 */
void FindPerm(vector<size_t> &perm, vector<size_t> &inv_perm, Graph &graph);

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
SPM *LoadMMF_RCM(const char *file_name, const RuntimeContext &rt_config);
void ReorderMatMMF(MMF &mat, const vector<size_t> &perm);
Graph& ConstructGraphMMF(Graph &graph, MMF &mat);

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
SPM *LoadCSR_RCM(IndexType *rowptr, IndexType *colind, ValueType *values,
                 size_t nr_rows, size_t nr_cols, bool zero_based,
                 const RuntimeContext &rt_config);
template<typename IndexType, typename ValueType>
void ReorderMatCSR(CSR<IndexType, ValueType> &mat, const vector<size_t> &perm,
                   vector<size_t> &inv_perm);
template<typename IterT>
Graph& ConstructGraphCSR(Graph &graph, IterT &iter, const IterT &iter_end,
                         size_t nr_nzeros, bool symmetric);

#endif // RCM_H__
