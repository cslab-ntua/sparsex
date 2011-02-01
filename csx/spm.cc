#include "spm.h"
#include "mmf.h"

#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <boost/function.hpp>
#include <boost/foreach.hpp>

#define FOREACH BOOST_FOREACH

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include "../../debug/debug_user.h"

extern "C" {
    #include "dynarray.h"
}

namespace bll = boost::lambda;
using namespace csx;

namespace csx {

std::ostream &operator<<(std::ostream &os, const DeltaRLE::StatsVal &stats)
{
    os << "nnz: " << stats.nnz;
    return os;
}

std::ostream &operator<<(std::ostream &os, const DeltaRLE &p)
{
    os << " (";
    p.PrintOn(os);
    os << " type:" << p.type << ") ";
    return os;
}

std::ostream &operator<<(std::ostream &out, CooElem p)
{
    out << "(" << std::setw(2) << p.y << "," << std::setw(2) << p.x << ")";
    return out;
}

std::ostream &operator<<(std::ostream &out, const SpmCooElem e)
{
    out << static_cast<CooElem>(e);
    if (e.pattern != NULL) {
        out << "->[" << *(e.pattern) << "]";
        out << " vals:{ ";
        for (int i=0; i<e.pattern->GetSize(); i++)
            out << e.vals[i] << " ";
        out << "}";
    } else {
        out << "v: " << e.val;
    }
    
    return out;
}

std::ostream &operator<<(std::ostream &out, const SpmRowElem &elem)
{
    out << "x:" << elem.x;
    if (elem.pattern) {
        out << *(elem.pattern);
        out << " vals:{ ";
        for (int i=0; i<elem.pattern->GetSize(); i++)
            out << elem.vals[i] << " ";
        out << "}";
    } else {
        out << " v:" << elem.val;
    }
    
    return out;
}

std::ostream &operator<<(std::ostream &out, SPM::PntIter pi)
{
    out << "<" << std::setw(2) << pi.row_idx << "," << std::setw(2)
        << pi.elm_idx << ">";
    return out;
}

#define BLOCK_ROW_MAP_NAME(r)   pnt_map_bR ## r
#define BLOCK_ROW_RMAP_NAME(r)  pnt_rmap_bR ## r
#define BLOCK_COL_MAP_NAME(c)   pnt_map_bC ## c
#define BLOCK_COL_RMAP_NAME(c)  pnt_rmap_bC ## c

#define DEFINE_BLOCK_ROW_MAP_FN(r)                                      \
    static inline void BLOCK_ROW_MAP_NAME(r)                            \
        (const CooElem &src, CooElem &dst)                              \
    {                                                                   \
        uint64_t src_x = src.x;                                         \
        uint64_t src_y = src.y;                                         \
                                                                        \
        dst.y = (src_y - 1) / r + 1;                                    \
        dst.x = (src_y - 1) % r + r * (src_x - 1) + 1;                  \
    }

#define DEFINE_BLOCK_ROW_RMAP_FN(r)                                     \
    static inline void BLOCK_ROW_RMAP_NAME(r)                           \
        (const CooElem &src, CooElem &dst)                              \
    {                                                                   \
        uint64_t src_x = src.x;                                         \
        uint64_t src_y = src.y;                                         \
                                                                        \
        dst.y = r * (src_y - 1) + (src_x - 1) % r + 1;                  \
        dst.x = (src_x - 1) / r + 1;                                    \
    }

#define DEFINE_BLOCK_COL_MAP_FN(c)                                      \
    static inline void BLOCK_COL_MAP_NAME(c)                            \
        (const CooElem &src, CooElem &dst)                              \
    {                                                                   \
        pnt_map_V(src, dst);                                            \
        BLOCK_ROW_MAP_NAME(c)(src, dst);                                \
    }

#define DEFINE_BLOCK_COL_RMAP_FN(c)                                     \
    static inline void BLOCK_COL_RMAP_NAME(c)                           \
        (const CooElem &src, CooElem &dst)                              \
    {                                                                   \
        BLOCK_ROW_RMAP_NAME(c)(src, dst);                               \
        pnt_rmap_V(src, dst);                                           \
    }

static inline void pnt_map_V(const CooElem &src, CooElem &dst)
{
    uint64_t src_x = src.x;
    uint64_t src_y = src.y;
    
    dst.x = src_y;
    dst.y = src_x;
}

static inline void pnt_rmap_V(const CooElem &src, CooElem &dst)
{
    pnt_map_V(src, dst);
}

static inline void pnt_map_D(const CooElem &src, CooElem &dst, uint64_t nrows)
{
    uint64_t src_x = src.x;
    uint64_t src_y = src.y;
    
    assert(nrows + src_x - src_y > 0);
    dst.y = nrows + src_x - src_y;
    dst.x = (src_x < src_y) ? src_x : src_y;
}

static inline void pnt_rmap_D(const CooElem &src, CooElem &dst, uint64_t nrows)
{
    uint64_t src_x = src.x;
    uint64_t src_y = src.y;
    
    if (src_y < nrows) {
        dst.x = src_x;
        dst.y = nrows + src_x - src_y;
    } else {
        dst.y = src_x;
        dst.x = src_y + src_x - nrows;
    }
}

static inline void pnt_map_rD(const CooElem &src, CooElem &dst, uint64_t ncols)
{
    uint64_t src_x = src.x;
    uint64_t src_y = src.y;
    uint64_t dst_y;
    
    dst.y = dst_y = src_x + src_y - 1;
    dst.x = (dst_y <= ncols) ? src_y : ncols + 1 - src_x;
}

static inline void pnt_rmap_rD(const CooElem &src, CooElem &dst, uint64_t ncols)
{
    uint64_t src_x = src.x;
    uint64_t src_y = src.y;
    
    if (src_y <= ncols) {
        dst.y = src_x;
        dst.x = src_y + 1 - src_x;
    } else {
        dst.y = src_x + src_y - ncols;
        dst.x = ncols + 1 - src_x;
    }
    
    dst.y = src_y - dst.x + 1;
}

DEFINE_BLOCK_ROW_MAP_FN(2)
DEFINE_BLOCK_ROW_MAP_FN(3)
DEFINE_BLOCK_ROW_MAP_FN(4)
DEFINE_BLOCK_ROW_MAP_FN(5)
DEFINE_BLOCK_ROW_MAP_FN(6)
DEFINE_BLOCK_ROW_MAP_FN(7)
DEFINE_BLOCK_ROW_MAP_FN(8)

DEFINE_BLOCK_ROW_RMAP_FN(2)
DEFINE_BLOCK_ROW_RMAP_FN(3)
DEFINE_BLOCK_ROW_RMAP_FN(4)
DEFINE_BLOCK_ROW_RMAP_FN(5)
DEFINE_BLOCK_ROW_RMAP_FN(6)
DEFINE_BLOCK_ROW_RMAP_FN(7)
DEFINE_BLOCK_ROW_RMAP_FN(8)

DEFINE_BLOCK_COL_MAP_FN(2)
DEFINE_BLOCK_COL_MAP_FN(3)
DEFINE_BLOCK_COL_MAP_FN(4)
DEFINE_BLOCK_COL_MAP_FN(5)
DEFINE_BLOCK_COL_MAP_FN(6)
DEFINE_BLOCK_COL_MAP_FN(7)
DEFINE_BLOCK_COL_MAP_FN(8)

DEFINE_BLOCK_COL_RMAP_FN(2)
DEFINE_BLOCK_COL_RMAP_FN(3)
DEFINE_BLOCK_COL_RMAP_FN(4)
DEFINE_BLOCK_COL_RMAP_FN(5)
DEFINE_BLOCK_COL_RMAP_FN(6)
DEFINE_BLOCK_COL_RMAP_FN(7)
DEFINE_BLOCK_COL_RMAP_FN(8)

void MakeRowElem(const CooElem &p, SpmRowElem *ret)
{
    ret->x = p.x;
    ret->val = p.val;
    ret->pattern = NULL;
}

void MakeRowElem(const SpmCooElem &p, SpmRowElem *ret)
{
    ret->x = p.x;
    ret->val = p.val;
    ret->pattern = (p.pattern == NULL) ? NULL : (p.pattern)->Clone();
}

void MakeRowElem(const SpmRowElem &p, SpmRowElem *ret)
{
    ret->x = p.x;
    ret->val = p.val;
    ret->pattern = (p.pattern == NULL) ? NULL : (p.pattern)->Clone();
}

}

SpmRowElem *SPM::RowBegin(uint64_t ridx)
{
    assert(ridx < this->rowptr_size - 1);
    return &this->elems[this->rowptr[ridx]];
}

SpmRowElem *SPM::RowEnd(uint64_t ridx)
{
    assert(ridx < this->rowptr_size - 1);
    return &this->elems[this->rowptr[ridx + 1]];
}

template <typename IterT>
uint64_t SPM::SetElems(IterT &pi, const IterT &pnts_end, uint64_t first_row,
                       unsigned long limit, uint64_t nr_elems, uint64_t nrows,
                       SPM::Builder *SpmBld)
{
    SpmRowElem *elem;
    uint64_t row_prev, row;
    
    row_prev = first_row;
    for (; pi != pnts_end; ++pi) {
        row = (*pi).y;
        if (row != row_prev) {
            assert(row > row_prev);
            if (limit && SpmBld->GetElemsCnt() >= limit)
                break;
            SpmBld->NewRow(row - row_prev);
            row_prev = row;
        }
        
        elem = SpmBld->AllocElem();
        MakeRowElem(*pi, elem);
    }
    
    return this->elems_size;
}

template <typename IterT>
uint64_t SPM::SetElems(IterT &pi, const IterT &pnts_end, uint64_t first_row,
                       unsigned long limit, uint64_t nr_elems, uint64_t nrows)
{
    SPM::Builder *SpmBld = new SPM::Builder(this, nr_elems, nrows);
    
    SetElems(pi, pnts_end, first_row, limit, nr_elems, nrows, SpmBld);
    SpmBld->Finalize();
    delete SpmBld;
    return this->elems_size;
}

SPM *SPM::LoadMMF_mt(const char *mmf_file, const long nr)
{
    SPM *ret;
    std::ifstream mmf;

    mmf.open(mmf_file);
    ret = LoadMMF_mt(mmf, nr);
    mmf.close();
    return ret;
}

SPM *SPM::LoadMMF_mt(std::istream &in, const long nr)
{
    MMF mmf(in);
    return LoadMMF_mt(mmf, nr);
}

SPM *SPM::LoadMMF_mt(MMF &mmf, const long nr)
{
    SPM *ret, *spm;
    long limit, cnt, row_start;
    MMF::iterator iter = mmf.begin();
    MMF::iterator iter_end = mmf.end();
    
    ret = new SPM[nr];
    row_start = limit = cnt = 0;
    for (long i=0; i<nr; i++) {
        spm = ret + i;
        limit = (mmf.nnz - cnt) / (nr - i);
        spm->nnz = spm->SetElems(iter, iter_end, row_start + 1, limit);
        spm->nrows = spm->rowptr_size - 1;
        spm->ncols = mmf.ncols;
        spm->row_start = row_start;
        row_start += spm->nrows;
        spm->type = HORIZONTAL;
        cnt += spm->nnz;
    }
    
    assert((uint64_t)cnt == mmf.nnz);    
    return ret;
}

SPM *SPM::LoadMMF(std::istream &in)
{
    return LoadMMF_mt(in, 1);
}

SPM *SPM::LoadMMF(const char *mmf_file)
{
    SPM *ret;
    std::ifstream mmf;
    
    mmf.open(mmf_file);
    ret = LoadMMF(mmf);
    mmf.close();
    return ret;
}

void SPM::Print(std::ostream &out)
{
    SPM::PntIter p, p_start, p_end;
    
    p_start = this->PointsBegin();
    p_end = this->PointsEnd();
    for (p = p_start; p != p_end; ++p)
        out << " " << (*p);
    out << std::endl;
}

void SPM::PrintElems(std::ostream &out)
{
    SPM::PntIter p, p_start, p_end;
    uint64_t y0;
    static int cnt = 1;
    
    y0 = this->row_start;
    p_start = this->PointsBegin();
    p_end = this->PointsEnd();
    for (p=p_start; p != p_end; ++p) {
        if ((*p).pattern == NULL) {
            out << std::setiosflags(std::ios::scientific) << row_start + (*p).y
                << " " << (*p).x << " " << (*p).val << " cnt:" << cnt++ << "\n";
            continue;
        }
        
        boost::function<void (CooElem &p)> xform_fn, rxform_fn;
        
        DeltaRLE *pat;
        CooElem start;
        double *vals;
        
        pat = (*p).pattern;
        start = static_cast<CooElem>(*p);
        vals = start.vals;
        std::cout << SpmTypesNames[pat->type] << std::endl;
        xform_fn = GetTransformFn(this->type, pat->type);
        rxform_fn = GetTransformFn(pat->type, this->type);
        if (xform_fn)
            xform_fn(start);
            
        DeltaRLE::Generator *g = (*p).pattern->generator(start);
        
        while (!g->IsEmpty()) {
            CooElem e = g->Next();
            if (rxform_fn)
                rxform_fn(e);
            out << std::setiosflags(std::ios::scientific) << row_start + e.y
                << " " << e.x << " " << *vals++ << " cnt:" << cnt++ << "\n";
        }
        
        out << "=== END OF PATTERN ===" << std::endl;
    }
}

void SPM::PrintStats(std::ostream& out)
{
    uint64_t nr_rows_with_patterns;
    uint64_t nr_patterns, nr_patterns_before;
    uint64_t nr_nzeros_block;
    uint64_t nr_transitions;
    uint64_t nr_xform_patterns[XFORM_MAX];

    nr_rows_with_patterns = this->GetNrRows();
    nr_patterns = 0;
    nr_nzeros_block = 0;
    nr_transitions = 0;
    memset(nr_xform_patterns, 0, sizeof(nr_xform_patterns));
    for (uint64_t i = 0; i < this->GetNrRows(); i++) {
        uint64_t pt_size, pt_size_before;
        SpmIterOrder pt_type, pt_type_before;
        
        nr_patterns_before = nr_patterns;
        
        const SpmRowElem *elem = this->RowBegin(i);
        
        if (elem->pattern) {
            pt_size_before = pt_size = elem->pattern->GetSize();
            pt_type_before = pt_type = elem->pattern->type;
        } else {
            pt_size_before = pt_size = 0;
            pt_type_before = pt_type = NONE;
        }

        for ( ; elem != this->RowEnd(i); elem++) {
            if (elem->pattern) {
                ++nr_patterns;
                pt_size = elem->pattern->GetSize();
                pt_type = elem->pattern->type;
                nr_nzeros_block += pt_size;
                if ((pt_type != pt_type_before) ||
                    (pt_size_before && pt_size != pt_size_before))
                    ++nr_transitions;
                ++nr_xform_patterns[elem->pattern->type];
                pt_size_before = pt_size;
            } else {
                pt_type = NONE;
                ++nr_xform_patterns[NONE];
                if (pt_type != pt_type_before)
                    ++nr_transitions;
            }

            pt_type_before = pt_type;
        }

        if (nr_patterns == nr_patterns_before)
            --nr_rows_with_patterns;
    }

    int nr_encoded_types = 0;
    
    for (int i = 1; i < XFORM_MAX; i++)
        if (nr_xform_patterns[i]) {
            ++nr_encoded_types;
            out << SpmTypesNames[i] << ": "
            << nr_xform_patterns[i] << std::endl;
        }

    out << "Encoded types = " << nr_encoded_types << ", " 
        << "Avg patterns/row = " 
        << nr_patterns / (double) nr_rows_with_patterns << ", " 
        << "Avg nonzeros/pattern = "
        << nr_nzeros_block / (double) nr_patterns << ", " 
        << " Avg pattern transitions/row = "
        << nr_transitions / (double) nr_rows_with_patterns
        << std::endl;
}

SPM::PntIter SPM::PointsBegin(uint64_t ridx)
{
    return PntIter(this, ridx);
}

SPM::PntIter SPM::PointsEnd(uint64_t ridx)
{
    if (ridx == 0)
        ridx = this->GetNrRows();
    return PntIter(this, ridx);
}

static inline bool elem_cmp_less(const SpmCooElem &e0, const SpmCooElem &e1)
{
    int ret;
    
    ret = CooCmp(static_cast<CooElem>(e0), static_cast<CooElem>(e1));
    return (ret < 0);
}

void SPM::Transform(SpmIterOrder t, uint64_t rs, uint64_t re)
{
    PntIter p0, pe, p;
    std::vector<SpmCooElem> elems;
    std::vector<SpmCooElem>::iterator e0, ee;
    boost::function<void (CooElem &p)> xform_fn;

    if (this->type == t)
        return;
    xform_fn = this->GetTransformFn(this->type, t);
    elems.reserve(this->elems_size);
    p0 = PointsBegin(rs);
    pe = PointsEnd(re);
    for(p=p0; p != pe; ++p) {
        SpmCooElem p_new = SpmCooElem(*p);
        xform_fn(p_new);
        elems.push_back(p_new);
    }
    
    e0 = elems.begin();
    ee = elems.end();
    sort(e0, ee, elem_cmp_less);
    SetElems(e0, ee, rs + 1);
    elems.clear();
    this->type = t;
}

inline TransformFn SPM::GetRevXformFn(SpmIterOrder type)
{
    boost::function<void (CooElem &p)> ret;
    
    switch(type) {
        case HORIZONTAL:
            break;
        case VERTICAL:
            ret = bll::bind(pnt_rmap_V, bll::_1, bll::_1);
            break;
        case DIAGONAL:
            ret = bll::bind(pnt_rmap_D, bll::_1, bll::_1, this->nrows);
            break;
        case REV_DIAGONAL:
            ret = bll::bind(pnt_rmap_rD, bll::_1, bll::_1, this->ncols);
            break;
        case BLOCK_ROW_TYPE_NAME(2):
            ret = bll::bind(BLOCK_ROW_RMAP_NAME(2), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(3):
            ret = bll::bind(BLOCK_ROW_RMAP_NAME(3), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(4):
            ret = bll::bind(BLOCK_ROW_RMAP_NAME(4), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(5):
            ret = bll::bind(BLOCK_ROW_RMAP_NAME(5), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(6):
            ret = bll::bind(BLOCK_ROW_RMAP_NAME(6), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(7):
            ret = bll::bind(BLOCK_ROW_RMAP_NAME(7), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(8):
            ret = bll::bind(BLOCK_ROW_RMAP_NAME(8), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(2):
            ret = bll::bind(BLOCK_COL_RMAP_NAME(2), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(3):
            ret = bll::bind(BLOCK_COL_RMAP_NAME(3), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(4):
            ret = bll::bind(BLOCK_COL_RMAP_NAME(4), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(5):
            ret = bll::bind(BLOCK_COL_RMAP_NAME(5), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(6):
            ret = bll::bind(BLOCK_COL_RMAP_NAME(6), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(7):
            ret = bll::bind(BLOCK_COL_RMAP_NAME(7), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(8):
            ret = bll::bind(BLOCK_COL_RMAP_NAME(8), bll::_1, bll::_1);
            break;
        default:
            std::cerr << "Unknown type: " << type << std::endl;
            assert(false);
    }
    
    return ret;
}

inline TransformFn SPM::GetXformFn(SpmIterOrder type)
{
    boost::function<void (CooElem &p)> ret;
    
    switch(type) {
        case HORIZONTAL:
            ret = NULL;
            break;
        case VERTICAL:
            ret = bll::bind(pnt_map_V, bll::_1, bll::_1);
            break;
        case DIAGONAL:
            ret = bll::bind(pnt_map_D, bll::_1, bll::_1, this->nrows);
            break;
        case REV_DIAGONAL:
            ret = bll::bind(pnt_map_rD, bll::_1, bll::_1, this->ncols);
            break;
        case BLOCK_ROW_TYPE_NAME(2):
            ret = bll::bind(BLOCK_ROW_MAP_NAME(2), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(3):
            ret = bll::bind(BLOCK_ROW_MAP_NAME(3), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(4):
            ret = bll::bind(BLOCK_ROW_MAP_NAME(4), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(5):
            ret = bll::bind(BLOCK_ROW_MAP_NAME(5), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(6):
            ret = bll::bind(BLOCK_ROW_MAP_NAME(6), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(7):
            ret = bll::bind(BLOCK_ROW_MAP_NAME(7), bll::_1, bll::_1);
            break;
        case BLOCK_ROW_TYPE_NAME(8):
            ret = bll::bind(BLOCK_ROW_MAP_NAME(8), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(2):
            ret = bll::bind(BLOCK_COL_MAP_NAME(2), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(3):
            ret = bll::bind(BLOCK_COL_MAP_NAME(3), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(4):
            ret = bll::bind(BLOCK_COL_MAP_NAME(4), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(5):
            ret = bll::bind(BLOCK_COL_MAP_NAME(5), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(6):
            ret = bll::bind(BLOCK_COL_MAP_NAME(6), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(7):
            ret = bll::bind(BLOCK_COL_MAP_NAME(7), bll::_1, bll::_1);
            break;
        case BLOCK_COL_TYPE_NAME(8):
            ret = bll::bind(BLOCK_COL_MAP_NAME(8), bll::_1, bll::_1);
            break;
        default:
            std::cerr << "Unknown type: " << type << std::endl;
            assert(false);
    }
    
    return ret;
}

inline TransformFn SPM::GetTransformFn(SpmIterOrder from_type,
                                       SpmIterOrder to_type)
{
    boost::function<void (CooElem &p)> xform_fn, rxform_fn;

    rxform_fn = GetRevXformFn(from_type);
    xform_fn = GetXformFn(to_type);
    if (xform_fn == NULL)
        return rxform_fn;
    if (rxform_fn != NULL)
        xform_fn = bll::bind(xform_fn,(bll::bind(rxform_fn, bll::_1), bll::_1));
    return xform_fn;
}

SPM *SPM::GetWindow(uint64_t rs, uint64_t length)
{
    if (rs + length > this->rowptr_size - 1)
        length = this->rowptr_size - rs - 1;

    SPM *ret = new SPM();
    uint64_t es = this->rowptr[rs];
    uint64_t ee = this->rowptr[rs+length];

    ret->rowptr = new uint64_t[length+1];
    ret->rowptr_size = length + 1;
    for (uint64_t i = 0; i < ret->rowptr_size; i++)
        ret->rowptr[i] = this->rowptr[rs+i] - es;
    ret->elems = &this->elems[es];
    ret->elems_size = ee - es;
    ret->nrows = ret->rowptr_size - 1;
    ret->ncols = this->ncols;
    ret->nnz = ret->elems_size;
    ret->row_start = this->row_start + rs;
    ret->type = this->type;
    ret->elems_mapped_ = true;
    assert(ret->rowptr[ret->rowptr_size-1] == ret->elems_size);
    return ret;
}

SPM *SPM::ExtractWindow(uint64_t rs, uint64_t length)
{
    std::vector<SpmCooElem> elems;
    std::vector<SpmCooElem>::iterator elem_begin, elem_end;
    PntIter p0, pe, p;
    SPM *ret = new SPM();
    
    if (rs + length > this->rowptr_size - 1)
        length = this->rowptr_size - rs - 1;
    
    uint64_t es = this->rowptr[rs];
    uint64_t ee = this->rowptr[rs+length];
    
    assert(es <= ee);  
    elems.reserve(ee - es);
    p0 = PointsBegin(rs);
    pe = PointsEnd(rs + length - 1);
    for (p = p0; p != pe; ++p)
        elems.push_back(*p);
    elem_begin = elems.begin();
    elem_end = elems.end();
    ret->SetElems(elem_begin, elem_end, rs+1);
    elems.clear();
    ret->nrows = ret->rowptr_size - 1;
    ret->ncols = this->ncols;
    ret->nnz = ret->elems_size;
    ret->row_start = this->row_start + rs;
    ret->type = this->type;
    ret->elems_mapped_ = true;   
    assert(ret->rowptr[ret->rowptr_size-1] == ret->elems_size);
    return ret;
}

void SPM::PutWindow(const SPM *window)
{
    assert(window);
    assert(this->type == window->type);

    uint64_t rs = window->row_start - this->row_start;
    uint64_t es = this->rowptr[rs];

    for (uint64_t i = 0; i < window->rowptr_size; i++)
        this->rowptr[rs+i] = es + window->rowptr[i];

    if (!window->elems_mapped_)
        memcpy(&this->elems[es], window->elems,
               window->elems_size*sizeof(SpmRowElem));
}

SpmRowElem *SPM::Builder::AllocElem()
{
    if (this->spm_->elems_mapped_)
        assert(dynarray_size(this->da_elems_) < this->spm_->elems_size);
    return (SpmRowElem *) dynarray_alloc(this->da_elems_);
}

SpmRowElem *SPM::Builder::AllocElems(uint64_t nr)
{
    if (this->spm_->elems_mapped_)
        assert(dynarray_size(this->da_elems_) + nr < this->spm_->elems_size);
    return (SpmRowElem *)dynarray_alloc_nr(this->da_elems_, nr);
}

uint64_t SPM::Builder::GetElemsCnt()
{
    return dynarray_size(this->da_elems_);
}

void SPM::Builder::NewRow(uint64_t rdiff)
{
    uint64_t *rowptr;
    uint64_t elems_cnt;
    
    elems_cnt = GetElemsCnt();
    rowptr = (uint64_t *)dynarray_alloc_nr(this->da_rowptr_, rdiff);
    for (uint64_t i=0; i < rdiff; i++)
        rowptr[i] = elems_cnt;
}

void SPM::Builder::Finalize()
{
    uint64_t *last_rowptr;
    
    last_rowptr = (uint64_t *)dynarray_get_last(this->da_rowptr_);
    if (*last_rowptr != dynarray_size(this->da_elems_))
        NewRow();
    if (!spm_->elems_mapped_ && spm_->elems)
	free(spm_->elems);
    if (spm_->rowptr)
	free(spm_->rowptr);
    if (spm_->elems_mapped_)
        assert(spm_->elems_size == dynarray_size(da_elems_));
	spm_->elems_size = dynarray_size(da_elems_);
    if (spm_->elems_mapped_)
        free(this->da_elems_);
    else
        spm_->elems = (SpmRowElem *)dynarray_destroy(da_elems_);
    
    spm_->rowptr_size = dynarray_size(da_rowptr_);
    spm_->rowptr = (uint64_t *)dynarray_destroy(da_rowptr_);
    this->da_elems_ = this->da_rowptr_ = NULL;
}

//*** Do I keep the code from now on;

void PrintTransform(uint64_t y, uint64_t x, TransformFn xform_fn,
                    std::ostream &out)
{
    uint64_t i,j;
    CooElem p;
    
    for (i=1; i <= y; i++) {
        for (j=1; j <= x;  j++) {
            p.y = i;
            p.x = j;
            xform_fn(p);
            out << p << " ";
        }
        
        out << std::endl;
    }
}

void PrintDiagTransform(uint64_t y, uint64_t x, std::ostream &out)
{
    boost::function<void (CooElem &p)> xform_fn;
    
    xform_fn = bll::bind(pnt_map_D, bll::_1, bll::_1, y);
    PrintTransform(y, x, xform_fn, out);
}

void PrintRDiagTransform(uint64_t y, uint64_t x, std::ostream &out)
{
    boost::function<void (CooElem &p)> xform_fn;
    
    xform_fn = bll::bind(pnt_map_rD, bll::_1, bll::_1, x);
    PrintTransform(y, x, xform_fn, out);
}

void TestTransform(uint64_t y, uint64_t x, TransformFn xform_fn, TransformFn rxform_fn)
{
	uint64_t i,j;
	CooElem p0, p1;
	for (i=1; i <= y; i++){
		for (j=1; j <= x;  j++){
			p0.y = i;
			p0.x = j;
			xform_fn(p0);
			p1.y = p0.y;
			p1.x = p0.x;
			rxform_fn(p1);
			if ( (p1.y != i) || (p1.x != j) ){
				std::cout << "Error for " << i << "," << j << std::endl;
				std::cout << "Transformed: " << p0.y << "," << p0.x << std::endl;
				std::cout << "rTransformed:" << p1.y << "," << p1.x << std::endl;
				exit(1);
			}
		}
	}
}

void TestRDiagTransform(uint64_t y, uint64_t x)
{
	boost::function<void (CooElem &p)> xform_fn, rxform_fn;
	xform_fn = bll::bind(pnt_map_rD, bll::_1, bll::_1, x);
	rxform_fn = bll::bind(pnt_rmap_rD, bll::_1, bll::_1, x);
	TestTransform(y, x, xform_fn, rxform_fn);
}


#if 0
int main(int argc, char **argv)
{
	PrintRDiagTransform(5, 5, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(10, 5, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(5, 10, std::cout);
	std::cout << std::endl;

	TestRDiagTransform(5,5);
	TestRDiagTransform(10,5);
	TestRDiagTransform(5,10);
}



void printXforms()
{
	std::cout << "Diagonal" << std::endl;
	PrintDiagTransform(10, 10, std::cout);
	std::cout << std::endl;
	PrintDiagTransform(5, 10, std::cout);
	std::cout << std::endl;
	PrintDiagTransform(10, 5, std::cout);
	std::cout << std::endl;
	std::cout << "Reverse Diagonal" << std::endl;
	PrintRDiagTransform(10, 10, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(5, 10, std::cout);
	std::cout << std::endl;
	PrintRDiagTransform(10, 5, std::cout);
	std::cout << std::endl;
}

void TestDiagTransform(uint64_t y, uint64_t x)
{
	boost::function<void (point_t &p)> xform_fn, rxform_fn;
	xform_fn = boost::bind(pnt_map_D, _1, _1, y);
	rxform_fn = boost::bind(pnt_rmap_D, _1, _1, y);
	TestTransform(y, x, xform_fn, rxform_fn);
}

void TestXforms()
{
	TestDiagTransform(10, 10);
	TestDiagTransform(5, 10);
	TestDiagTransform(10, 5);
	TestRDiagTransform(10, 10);
	TestRDiagTransform(5, 10);
	TestRDiagTransform(10, 5);
}
#endif

#if 0
void test_drle()
{
	std::vector<int> v_in, delta;
	std::vector< RLE<int> > drle;

	v_in.resize(16);
	for (int i=0; i<16; i++){
		v_in[i] = i;
	}
	FOREACH(int &v, v_in){
		std::cout << v << " ";
	}
	std::cout << std::endl;

	delta = DeltaEncode(v_in);
	FOREACH(int &v, delta){
		std::cout << v << " ";
	}
	std::cout << std::endl;

	drle = RLEncode(delta);RMAP_
	FOREACH(RLE<int> &v, drle){
		std::cout << "(" << v.val << "," << v.freq << ") ";
	}
	std::cout << std::endl;
}
#endif

#if 0
int main(int argc, char **argv)
{
	SPM obj;
	obj.LoadMMF();
	//obj.DRLEncode();
	//obj.Print();
	//std::cout << obj.rows;
	//obj.Print();
	//obj.Draw("1.png");
	//std::cout << "\n";
	//obj.Transform(VERTICAL);
	//std::cout << obj.rows;
	//obj.Print();
	//obj.DRLEncode();
	//std::cout << obj.rows;
}
#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
