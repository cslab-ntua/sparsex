#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "reduce_bandwidth.h"

using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS, 
                       property<vertex_color_t, default_color_type,
                                property<vertex_degree_t,int> > > Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertices_size_type size_type;
typedef std::pair<std::size_t, std::size_t> Pair;

Matrix::Matrix(uint64_t n, uint64_t nnz)
{
    n_ = n;
    nnz_ = nnz;
    elems_.reserve(nnz);
}

void Matrix::Insert(uint64_t row, uint64_t col, double val)
{
    elem_t elem;
    
    elem.row = row;
    elem.col = col;
    elem.val = val;
    elems_.push_back(elem);
}

void Matrix::Get(uint64_t index, uint64_t *row, uint64_t *col, double *val)
{
    *row = elems_[index].row;
    *col = elems_[index].col;
    *val = elems_[index].val;
}

void Matrix::Update(uint64_t index, uint64_t row, uint64_t col)
{
    elems_[index].row = row;
    elems_[index].col = col;
}

void Matrix::Sort()
{
    sort(elems_.begin(), elems_.end(), Compare);
}

void Matrix::Print()
{   
    cout << "Elements of Matrix" << endl;
    cout << "------------------" << endl;
    for (uint64_t i = 0; i < nnz_; i++)
        cout << elems_[i].row << " " << elems_[i].col << " " << elems_[i].val <<
                endl;
}

int main(int argc, char **argv)
{
    uint64_t nrows, ncols, nnz, nr_edges;
    uint64_t row, col;
    double val;
    Matrix *matrix;
    uint64_t ob, nb;

    uint64_t i;
    char line[512];

    if (argc != 3) {
        cerr << "Wrong Input" << endl << "Correct Format: ./main <old_filename>"
             << " <new_filename>" << endl;
        return 0;   
    }

    ifstream fin (argv[1]);
    assert(fin.is_open() && "Unable to open input file\n");
    
    fin.getline(line, sizeof(line));
    sscanf(line, "%lu %lu %lu", &nrows, &ncols, &nnz);
    matrix = new Matrix(nrows, nnz);
    nr_edges = (nnz - nrows) / 2;

    fin.getline(line, sizeof(line));
    while (!fin.eof()) {
        sscanf(line, "%lu %lu %lf", &row, &col, &val);
        matrix->Insert(row, col, val);
        fin.getline(line, sizeof(line));
    }
    assert(row == nrows);
    fin.close();
    
    // matrix->Print();
    
    Pair *edges = (Pair *) malloc(nr_edges*sizeof(Pair));
    
    uint64_t index = 0;
    for (i = 0; i < nnz; i++) {
        matrix->Get(i, &row, &col, &val);

        if (row < col)
            edges[index++] = Pair(row-1, col-1);
    }
    assert(index == nr_edges);

    Graph G(nrows);
    for (i = 0; i < nr_edges; ++i)
        add_edge(edges[i].first, edges[i].second, G);
    graph_traits<Graph>::vertex_iterator ui, ui_end;
    property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
        deg[*ui] = degree(*ui, G);
        
    property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, G);
    
    ob = bandwidth(G);
    cout << "Original Bandwidth: " << ob << endl;
    vector<Vertex> inv_perm(num_vertices(G));
    vector<size_type> perm(num_vertices(G));

    //reverse cuthill_mckee_ordering
    cuthill_mckee_ordering(G, inv_perm.rbegin(), get(vertex_color, G),
                           make_degree_map(G));

    uint64_t *map = (uint64_t *) malloc(nnz*sizeof(uint64_t));
    
    i = 0;
    for (std::vector<Vertex>::const_iterator k=inv_perm.begin();
        k != inv_perm.end(); ++k) {
        map[index_map[*k]] = i;
        i++;
    }

    /*cout << "Reverse Cuthill-McKee ordering: ";
    for (i = 0; i < nrows; ++i)
        cout << map[i] << " ";
    cout << endl;*/

    for (size_type c = 0; c != inv_perm.size(); ++c)
        perm[index_map[inv_perm[c]]] = c;
    
    nb = bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0]));
    cout << "Final Bandwidth: " << nb << endl;

    for (i = 0; i < nnz; i++) {
        matrix->Get(i, &row, &col, &val);
        matrix->Update(i, map[row-1]+1, map[col-1]+1);
    }

    matrix->Sort();

    ofstream fout (argv[2]);
    assert(fout.is_open() && "Unable to open output file\n");
    fout << nrows << " " << ncols << " " << nnz << endl;
    for (i = 0; i < nnz; i++) {
        matrix->Get(i, &row, &col, &val);
        fout << row << " " << col << " " << val << endl;
    }
    fout.close();
    
    // matrix->Print();

    free(map);
    delete matrix;

    return 0;
}
