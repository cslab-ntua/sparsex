#ifndef RCM_IMPL_H__
#define RCM_IMPL_H__

//#include "rcm.h"
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
    perm.reserve(num_vertices(graph));
    invperm.reserve(num_vertices(graph));

    // Reverse Cuthill Mckee Ordering
    cuthill_mckee_ordering(graph, inv_perm.rbegin(), get(vertex_color, graph),
                           make_degree_map(graph));

    for (size_t i = 0; i != inv_perm.size(); ++i) {
        perm[index_map[inv_perm[i]]] = i;
        invperm.push_back(inv_perm[i]);
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

template<typename IterT>
Graph& ConstructGraphCSR(Graph &graph, IterT &iter, const IterT &iter_end,
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
void ReorderMatCSR(CSR<IndexType, ValueType> &mat, const vector<size_t> &perm,
                   vector<size_t> &inv_perm)
{
    // Apply permutation only to colind
    for (size_t i = 0; i < mat.GetNrNonzeros(); i++) {
        mat.colind_[i] = perm[mat.colind_[i] - !mat.IsZeroBased()] +
            !mat.IsZeroBased();
    }

    // Sort colind and values per row
    size_t row_start, length;
    for (size_t i = 0; i < mat.GetNrRows(); i++) {
        row_start = mat.rowptr_[i] - !mat.IsZeroBased();
        length = mat.rowptr_[i+1] - mat.rowptr_[i];
        std::sort(get_CSR_iterator<IndexType*, ValueType*>
                  (&mat.colind_[row_start], &mat.values_[row_start]), 
                  get_CSR_iterator<IndexType*, ValueType*>
                  (&mat.colind_[row_start] + length, &mat.values_[row_start] + length), 
                  CSR_Comp<IndexType*, ValueType*>());
    }
    mat.SetReordered(inv_perm);
    //assert(inv_perm.capacity() == 0);
}

template<typename IndexType, typename ValueType>
SPM *LoadCSR_RCM(IndexType *rowptr, IndexType *colind, ValueType *values,
                 size_t nr_rows, size_t nr_cols, bool zero_based,
                 const RuntimeContext &rt_config)
{
    SPM *spms = NULL;
    CSR<IndexType, ValueType> mat(rowptr, colind, values, 
                                  nr_rows, nr_cols, zero_based);

    typename CSR<IndexType, ValueType>::iterator iter = mat.begin();
    typename CSR<IndexType, ValueType>::iterator iter_end = mat.end();
    vector<size_t> perm, inv_perm;

    cout << "=== START REORDERING ===\n";

    // Construct graph
    Graph graph(nr_rows);
    graph = ConstructGraphCSR(graph, iter, iter_end, mat.GetNrNonzeros(),
                              mat.IsSymmetric());
    // Find permutation
    FindPerm(perm, inv_perm, graph);
    // Reorder original matrix
    ReorderMatCSR<IndexType, ValueType>(mat, perm, inv_perm);
    // Create internal representation
    spms = SPM::DoLoadMatrix<CSR<IndexType, ValueType> >
        (mat, rt_config.GetNrThreads());
    //spms->Print();

    cout << "=== END REORDERING ===\n";

    return spms;
}

#endif // RCM_IMPL_H__
