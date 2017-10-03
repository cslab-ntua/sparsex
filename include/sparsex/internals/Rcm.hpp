/*
 * Copyright (C) 2012-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Rcm.hpp
 * \brief Reverse Cuthill-Mckee Ordering Interface
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_RCM_HPP
#define SPARSEX_INTERNALS_RCM_HPP

#include <sparsex/internals/Csr.hpp>
#include <sparsex/internals/Mmf.hpp>
#include <sparsex/internals/logger/Logger.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

using namespace sparsex::io;

namespace sparsex {
  namespace utilities {

    using namespace boost;

#define bad_reorder 0

    typedef pair<size_t, size_t> Pair;
    typedef adjacency_list<vecS, vecS, undirectedS,
			   property<vertex_color_t, default_color_type,
				    property<vertex_degree_t,int> > > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::vertices_size_type size_type;

    template<typename IndexType, typename ValueType>
    struct ColumnCompare {
    public:
      bool operator() (const pair<IndexType, ValueType> &lhs,
		       const pair<IndexType, ValueType> &rhs) const
      {
	if (lhs.first < rhs.first) return true;
	if (lhs.first > rhs.first) return false;
	return false;
      }
    };

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
     *  @return            spm class object with the characteristics of the
     *                     matrix.
     */
    template<typename IndexType, typename ValueType>
    void ReorderMat_MMF(MMF<IndexType, ValueType>& mat,
			const vector<size_t>& perm);
    template<typename IndexType, typename ValueType>
    Graph& ConstructGraph_MMF(Graph& graph, MMF<IndexType, ValueType>& mat);
    template<typename IndexType, typename ValueType>
    void DoReorder_RCM(MMF<IndexType, ValueType>& mat, vector<size_t>& perm);

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
     *  @return            spm class object with the characteristics of the
     *                     matrix.
     */
    template<typename IndexType, typename ValueType>
    void ReorderMat_CSR(CSR<IndexType, ValueType>& mat,
			const vector<size_t>& perm,
			vector<size_t>& inv_perm);
    template<typename IterT>
    Graph& ConstructGraph_CSR(Graph& graph, IterT& iter, const IterT& iter_end,
			      size_t nr_nzeros, bool symmetric);
    template<typename IndexType, typename ValueType>
    void DoReorder_RCM(CSR<IndexType, ValueType>& mat, vector<size_t>& perm);

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
      LOG_INFO << "Original Bandwidth: " << ob << "\n";

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
      LOG_INFO << "Final Bandwidth: " << nb << "\n";
    }

    template<typename IndexType, typename ValueType>
    Graph& ConstructGraph_MMF(Graph& graph, MMF<IndexType, ValueType>& mat)
    {
      // The flag must be set before using the MMF iterator
      mat.SetReordered(true);

      // Make a pessimistic guess for nr_edges
      size_t nr_edges;
      if (mat.IsSymmetric()) {
        // If main diagonal is full (nr_nzeros - nr_rows) / 2
        nr_edges = mat.GetNrNonzeros() / 2;
      } else {
        nr_edges = mat.GetNrNonzeros();
      }

      Pair *edges = new Pair[nr_edges];

      typename MMF<IndexType, ValueType>::iterator iter = mat.begin();
      typename MMF<IndexType, ValueType>::iterator iter_end = mat.end();
      size_t index = 0;
      if (!mat.IsSymmetric() && !mat.IsColWise()) {
        mat.InitMatrix(mat.GetNrNonzeros());
        for (; iter != iter_end; ++iter) {
	  mat.InsertElement(*iter);
	  if ((*iter).GetRow() != (*iter).GetCol())
	    edges[index++] = make_pair((*iter).GetRow() - 1,
				       (*iter).GetCol() - 1);
        }
      } else {
        for (;iter != iter_end; ++iter)
	  if ((*iter).GetRow() < (*iter).GetCol())
	    edges[index++] = make_pair((*iter).GetRow() - 1,
				       (*iter).GetCol() - 1);
      }

      // index -> actual nr_edges
      assert(index <= nr_edges);
      if (index == 0) {
        LOG_WARNING << "no reordering available for this matrix\n";
        delete[] edges;
        throw bad_reorder;
      }

      for (size_t i = 0; i < index; ++i) {
        add_edge(edges[i].first, edges[i].second, graph);
      }

      delete[] edges;
      return graph;
    }

    template<typename IndexType, typename ValueType>
    void ReorderMat_MMF(MMF<IndexType, ValueType>& mat,
			const vector<size_t>& perm)
    {
      for (size_t i = 0; i < mat.GetNrNonzeros(); i++) {
        pair<IndexType, ValueType> coord = mat.GetCoordinates(i);
        mat.SetCoordinates(i, perm[coord.first-1] + 1,
                           perm[coord.second-1] + 1);
      }

      mat.Sort();
    }

    template<typename IndexType, typename ValueType>
    void DoReorder_RCM(MMF<IndexType, ValueType>& mat, vector<size_t> &perm)
    {
      vector<size_t> inv_perm;

      LOG_INFO << "Reordering input matrix...\n";
      // Construct graph
      Graph graph(mat.GetNrRows());
      try {
        graph = ConstructGraph_MMF(graph, mat);
      } catch (int e) {
        mat.ResetStream();
        mat.SetReordered(false);
        LOG_WARNING << "reordering failed\n";
        return;
      }
      // Find permutation
      FindPerm(perm, inv_perm, graph);
      // Reorder original matrix
      ReorderMat_MMF(mat, perm);
      LOG_INFO << "Reordering complete\n";
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

      Pair *edges = new Pair[nr_edges];
      size_t index = 0;
      if (symmetric) {
        for (;iter != iter_end; ++iter) {
	  if ((*iter).GetRow() < (*iter).GetCol()) {
	    edges[index++] = make_pair(((*iter).GetRow()) - 1,
				       ((*iter).GetCol()) - 1);
	  }
        }
      } else {
        for (;iter != iter_end; ++iter) {
	  if ((*iter).GetRow() != (*iter).GetCol()) {
	    edges[index++] = make_pair(((*iter).GetRow()) - 1,
				       ((*iter).GetCol()) - 1);
	  }
        }
      }

      // index -> actual nr_edges
      assert(index <= nr_edges);

      if (index == 0) {
        LOG_WARNING << "no reordering available for this matrix\n";
        delete[] edges;
        throw bad_reorder;
      }

      for (size_t i = 0; i < index; i++) {
        add_edge(edges[i].first, edges[i].second, graph);
      }

      delete[] edges;
      return graph;
    }

    template<typename IndexType, typename ValueType>
    void ReorderMat_CSR(CSR<IndexType, ValueType>& mat,
			const vector<size_t>& perm,
			vector<size_t>& inv_perm)
    {
      // Have to work on a copy of colind and values if CSR structures are not
      // be shared
      IndexType *new_colind = new IndexType[mat.GetNrNonzeros()];
      ValueType *new_values = new ValueType[mat.GetNrNonzeros()];
      memcpy(new_values, mat.values_, sizeof(ValueType) * mat.GetNrNonzeros());

      // Apply permutation only to colind
      for (size_t i = 0; i < mat.GetNrNonzeros(); i++) {
        new_colind[i] = perm[mat.colind_[i] - !mat.IsZeroBased()] +
	  !mat.IsZeroBased();
      }

      mat.colind_ = new_colind;
      mat.values_ = new_values;

      // Simultaneously sort colind and values per row
      for (size_t i = 0; i < mat.GetNrRows(); i++) {
        sort(mat.row_begin(i), mat.row_end(i),
             ColumnCompare<IndexType, ValueType>());
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

      LOG_INFO << "Reordering input matrix...\n";
      // Construct graph
      Graph graph(mat.GetNrRows());
      try {
        graph = ConstructGraph_CSR(graph, iter, iter_end, mat.GetNrNonzeros(),
                                   mat.IsSymmetric());
      } catch (int e) {
        LOG_INFO << "reordering failed\n";
        return;
      }
      // Find permutation
      FindPerm(perm, inv_perm, graph);
      // Reorder original matrix
      ReorderMat_CSR<IndexType, ValueType>(mat, perm, inv_perm);
      LOG_INFO << "Reordering complete\n";
    }

  } // end of namespace utilities
} // end of namespace sparsex

#endif // SPARSEX_INTERNALS_RCM_HPP
