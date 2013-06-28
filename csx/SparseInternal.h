/* -*- C++ -*-
 *
 * SparseInternal.h --  Internal representation of sparse matrices.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SPARSEINTERNAL_H__
#define CSX_SPARSEINTERNAL_H__

#include "SparseUtil.h"
#include "SparsePartition.h"
#include "dynarray.h"

namespace csx {

template<typename IndexType, typename ValueType>
class SparseInternal
{
public:

    SparseInternal(IndexType nr_rows, IndexType nr_cols, IndexType nr_nzeros,
                   size_t nr_partitions)
        : nr_rows_(nr_rows),
          nr_cols_(nr_cols),
          nr_nzeros_(nr_nzeros),
          nr_partitions_(nr_partitions),
          partitions_(NULL)
    {}

    ~SparseInternal()
    {
        if (partitions_)
            delete[] partitions_;
    }

    IndexType GetNrRows() const
    {
        return nr_rows_;
    }

    IndexType GetNrCols() const
    {
        return nr_cols_;
    }

    IndexType GetNrNonzeros() const
    {
        return nr_nzeros_;
    }

    SparsePartition<IndexType, ValueType> *GetPartition(int id)
    {
        return partitions_ + id;
    }

    template<typename MatrixType>
    void BuildPool(MatrixType &mat)
    {
        typename MatrixType::iterator iter = mat.begin();
        typename MatrixType::iterator iter_end = mat.end();

        data_pool_.SetPool(iter, iter_end, mat.GetNrNonzeros());
    }

    void BuildPartitions(IndexType nr_nzeros, size_t nr);

    template<typename IterT>
    SparsePartition<IndexType, ValueType> *BuildPartition(IterT &iter,
                                                          const IterT &iter_end,
                                                          size_t id,
                                                          IndexType row_start,
                                                          IndexType elem_start,
                                                          unsigned long limit,
                                                          IndexType nr_rows,
                                                          IndexType nr_cols);

    template<typename MatrixType>
    static SparseInternal *DoLoadMatrix(MatrixType &mat, size_t nr)
    {
        SparseInternal *ret = new SparseInternal(mat.GetNrRows(),
                                                 mat.GetNrCols(),
                                                 mat.GetNrNonzeros(),
                                                 nr);
        ret->nr_rows_ = mat.GetNrRows();
        ret->nr_cols_ = mat.GetNrCols();
        ret->nr_nzeros_ = mat.GetNrNonzeros();
        ret->BuildPool(mat);
        ret->BuildPartitions(mat.GetNrNonzeros(), nr);
        return ret;
    }

    void Print(std::ostream &os) const
    {
        os << "Elements of Encoded Matrix" << endl;
        os << "--------------------------" << endl;
        for (size_t i = 0; i < nr_partitions_; i++) {
            os << "Partition " << i << ":" << endl;
            partitions_[i].Print(os);
            os << endl;
        }
    }

private:

    class ElemPool
    {
    public:

        ElemPool() : nr_elems_(0), data_(NULL) {}
        ElemPool(IndexType nr_elems) : nr_elems_(nr_elems), data_(NULL) {}

        ~ElemPool()
        {
            if (data_)
                free(data_);
        }

        IndexType GetNrNonzeros() const
        {
            return nr_elems_;
        }

        template<typename IterT>
        void SetPool(IterT& iter, const IterT& iter_end, IndexType nr_elems);

        Elem<IndexType, ValueType>& operator[](const IndexType idx)
        {
            // if (idx < 0 || idx >= nr_elems_) {
            //     std::cerr << "Index out of bounds!" << std::endl;
            //     exit(1);
            // }
            assert(idx >= 0 && idx < nr_nzeros_ && "Index out of bounds!");
            return data_[idx];
        }

        void PrintOn(std::ostream &out);

        class iterator;
        friend class iterator;
    
        iterator begin()
        {
            return iterator(this, 0);
        }

        iterator end()
        {
            return iterator(this, this->nr_elems_);
        }
 
    private:
        class Builder;
    public:
        IndexType nr_elems_;
        Elem<IndexType, ValueType> *data_;
    };

private:
    IndexType nr_rows_, nr_cols_, nr_nzeros_, nr_partitions_; 
    ElemPool data_pool_;
    SparsePartition<IndexType, ValueType> *partitions_;
};

/* Builder class */
template<typename IndexType, typename ValueType> 
class SparseInternal<IndexType, ValueType>::ElemPool::Builder
{
public:
    Builder(ElemPool *pool, size_t nr_elems = 0);
    ~Builder();
    Elem<IndexType, ValueType> *AllocElem();
    void Finalize();

private:
    ElemPool *pool_;
    dynarray_t *da_elems_;
};

/* Iterator class */
template<typename IndexType, typename ValueType>
class SparseInternal<IndexType, ValueType>::ElemPool::iterator
    : public std::iterator<std::forward_iterator_tag, Elem<IndexType,
                                                           ValueType> >
{
public:
    iterator() {}
    
    iterator(ElemPool *pool, IndexType cnt): pool_(pool), cnt_(cnt) {}
    
    bool operator==(const iterator &i)
    {
        return (pool_ == i.pool_) && (cnt_ == i.cnt_);
    }

    bool operator!=(const iterator &i)
    {
        return !(*this == i);
    }

    void operator++()
    {
        ++cnt_;
    }

    Elem<IndexType, ValueType> operator*()
    {
        return pool_->data_[cnt_];
    }

private:
    ElemPool *pool_;
    IndexType cnt_;
};

}  //csx namespace end

/*
 * SparseInternal class implementation
 */
template<typename IndexType, typename ValueType>
void SparseInternal<IndexType, ValueType>::BuildPartitions(IndexType nr_nzeros,
                                                           size_t nr)
{
    SparsePartition<IndexType, ValueType> *partition;
    unsigned long limit;
    uint64_t row_start, cnt;

    typename ElemPool::iterator iter = data_pool_.begin();
    typename ElemPool::iterator iter_end = data_pool_.end();

    partitions_ = new SparsePartition<IndexType, ValueType>[nr];
    row_start = limit = cnt = 0;
    for (size_t i = 0; i < nr; ++i) {
        limit = (nr_nzeros - cnt) / (nr - i);
        partition = BuildPartition(iter, iter_end, i, row_start, cnt,
                                   limit, nr_rows_, nr_cols_);
        row_start += partition->GetNrRows();
        cnt += partition->GetNrNonzeros();
    }
    assert((uint64_t) cnt == nr_nzeros_);
}

template<typename IndexType, typename ValueType>
template<typename IterT>
SparsePartition<IndexType, ValueType> *SparseInternal<IndexType, ValueType>::
BuildPartition(IterT &iter, const IterT &iter_end, size_t id,
               IndexType row_start, IndexType elem_start, unsigned long limit,
               IndexType nr_rows, IndexType nr_cols)
{
    SparsePartition<IndexType, ValueType> *ret;

    ret = partitions_ + id;
    IndexType nr_nzeros = ret->SetRows(iter, iter_end, row_start + 1,
                                        limit, nr_rows + 1);

    ret->SetNrNonzeros(nr_nzeros);
    ret->SetNrRows(ret->GetRowPtrSize() - 1);
    //ret->SetNrRows(nr_rows);
    //std::cout << ret->GetNrRows() << std::endl;
    ret->SetNrCols(nr_cols);
    ret->SetRowStart(row_start);
    ret->SetType(HORIZONTAL);
    ret->SetElemsMapped(true);
    ret->SetElems(data_pool_.data_ + elem_start);
    return ret;
}

template<typename IndexType, typename ValueType> 
template<typename IterT> 
void SparseInternal<IndexType, ValueType>::
ElemPool::SetPool(IterT& iter, const IterT& iter_end, IndexType nr_elems)
{
    Elem<IndexType, ValueType> *elem;
    Builder *Bld = new Builder(this, nr_elems);

    for (; iter != iter_end; ++iter) {
        elem = Bld->AllocElem();
        MakeRowElem<IndexType, ValueType>(*iter, elem);
    }
    Bld->Finalize();

    delete Bld;
    nr_elems_ = nr_elems;
}

template<typename IndexType, typename ValueType>
void SparseInternal<IndexType, ValueType>::ElemPool::PrintOn(ostream& os)
{
    for (size_t i = 0; i < nr_elems_; i++)
        os << data_[i] << " ";
    os << endl;
}

/*
 * Builder class implementation
 */
template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType>::ElemPool::Builder::Builder(ElemPool *pool,
                                                                 size_t nr_elems)
    : pool_(pool)
{
    da_elems_ = dynarray_create(sizeof(Elem<IndexType, ValueType>),
                                512, nr_elems);
}

template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType>::ElemPool::Builder::~Builder()
{
    assert(da_elems_ == NULL && "da_elems_ not destroyed");
}

template<typename IndexType, typename ValueType>
Elem<IndexType, ValueType> *SparseInternal<IndexType, ValueType>::
ElemPool::Builder::AllocElem()
{
    return (Elem<IndexType, ValueType> *) dynarray_alloc(da_elems_);
}

template<typename IndexType, typename ValueType>
void SparseInternal<IndexType, ValueType>::ElemPool::Builder::Finalize()
{
    pool_->data_ = (Elem<IndexType, ValueType> *) dynarray_destroy(da_elems_);
    da_elems_ = NULL;
}

#endif // CSX_SPARSEINTERNAL_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
