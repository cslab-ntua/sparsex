#ifndef REDUCE_BANDWIDTH_H_
#define REDUCE_BANDWIDTH_H_

extern "C" {
#include <inttypes.h>
}
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

typedef struct {
    uint64_t row;
    uint64_t col;
    double val;
} elem_t;

bool Compare(elem_t i, elem_t j)
{
    if (i.row < j.row)
        return true;
    if (i.row > j.row)
        return false;
    if (i.col < j.col)
        return true;
    return false;
}

class Matrix
{
private:
    uint64_t n_;
    uint64_t nnz_;
    std::vector<elem_t> elems_;
public:
    Matrix(uint64_t n, uint64_t nnz);
    ~Matrix() {};

    void Insert(uint64_t row, uint64_t col, double val);
    void Get(uint64_t index, uint64_t *row, uint64_t *col, double *val);
    void Update(uint64_t index, uint64_t row, uint64_t col);
    void Sort();
    
    void Print();
};

#endif
