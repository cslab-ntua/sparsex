/*
 * spm.cc --  Internal representation of sparse matrices (implementation)
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include "spm.h"

#define FOREACH BOOST_FOREACH

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
    os << " type:" << p.GetType() << ") ";
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
        for (int i = 0; i < e.pattern->GetSize(); i++)
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
    out << " in_pattern:" << elem.in_pattern;
    out << " pattern_start:" << elem.pattern_start;
    if (elem.pattern) {
        out << *(elem.pattern);
        out << " vals:{ ";
        for (int i = 0; i < elem.pattern->GetSize(); i++)
            out << elem.vals[i] << " ";
            
        out << "}";
    } else {
        out << " v:" << elem.val;
    }

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

DEFINE_BLOCK_ROW_MAP_FN(1)
DEFINE_BLOCK_ROW_MAP_FN(2)
DEFINE_BLOCK_ROW_MAP_FN(3)
DEFINE_BLOCK_ROW_MAP_FN(4)
DEFINE_BLOCK_ROW_MAP_FN(5)
DEFINE_BLOCK_ROW_MAP_FN(6)
DEFINE_BLOCK_ROW_MAP_FN(7)
DEFINE_BLOCK_ROW_MAP_FN(8)

DEFINE_BLOCK_ROW_RMAP_FN(1)
DEFINE_BLOCK_ROW_RMAP_FN(2)
DEFINE_BLOCK_ROW_RMAP_FN(3)
DEFINE_BLOCK_ROW_RMAP_FN(4)
DEFINE_BLOCK_ROW_RMAP_FN(5)
DEFINE_BLOCK_ROW_RMAP_FN(6)
DEFINE_BLOCK_ROW_RMAP_FN(7)
DEFINE_BLOCK_ROW_RMAP_FN(8)

DEFINE_BLOCK_COL_MAP_FN(1)
DEFINE_BLOCK_COL_MAP_FN(2)
DEFINE_BLOCK_COL_MAP_FN(3)
DEFINE_BLOCK_COL_MAP_FN(4)
DEFINE_BLOCK_COL_MAP_FN(5)
DEFINE_BLOCK_COL_MAP_FN(6)
DEFINE_BLOCK_COL_MAP_FN(7)
DEFINE_BLOCK_COL_MAP_FN(8)

DEFINE_BLOCK_COL_RMAP_FN(1)
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
    ret->in_pattern = p.in_pattern;
    ret->pattern_start = p.pattern_start;
    ret->pattern = (p.pattern == NULL) ? NULL : (p.pattern)->Clone();
}

void MakeRowElem(const SpmRowElem &p, SpmRowElem *ret)
{
    ret->x = p.x;
    ret->val = p.val;
    ret->in_pattern = p.in_pattern;
    ret->pattern_start = p.pattern_start;
    ret->pattern = (p.pattern == NULL) ? NULL : (p.pattern)->Clone();
}

}

SpmRowElem *SPM::RowBegin(uint64_t ridx)
{
    assert(ridx < rowptr_size_ - 1 && "ridx out of bounds");
    return &elems_[rowptr_[ridx]];
}

SpmRowElem *SPM::RowEnd(uint64_t ridx)
{
    assert(ridx < rowptr_size_ - 1 && "ridx out of bounds");
    return &elems_[rowptr_[ridx+1]];
}

SPM *SPM::LoadMMF_mt(const char *mmf_file, const long nr)
{
    SPM *ret;
    std::ifstream mmf;

    mmf.open(mmf_file);

    if (mmf.good()) {
        ret = LoadMMF_mt(mmf, nr);
    } else {
        std::cerr << "File error!" << std::endl;
        exit(1);
    }

    mmf.close();
    return ret;
}

/*SPM *SPM::LoadMMF_mt(std::istream &in, const long nr)
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
    for (long i = 0; i < nr; ++i) {
        spm = ret + i;
        limit = (mmf.GetNrNonzeros() - cnt) / (nr - i);
        spm->nr_nzeros_ = spm->SetElems(iter, iter_end, row_start + 1, limit,
                                        limit + mmf.GetNrRows() - 1,
                                        mmf.GetNrRows() + 1);
        spm->nr_rows_ = spm->rowptr_size_ - 1;
        spm->nr_cols_ = mmf.GetNrCols();
        spm->row_start_ = row_start;
        row_start += spm->nr_rows_;
        spm->type_ = HORIZONTAL;
        cnt += spm->nr_nzeros_;
    }

    assert((uint64_t)cnt == mmf.GetNrNonzeros());
    return ret;
}*/

SPM *SPM::LoadMMF(std::istream &in)
{
    return LoadMMF_mt(in, 1);
}

SPM *SPM::LoadMMF(const char *mmf_file)
{
    SPM *ret;
    std::ifstream mmf;

    mmf.open(mmf_file);

    if (mmf.good()) {
        ret = LoadMMF(mmf);
    } else {
        std::cerr << "File error!" << std::endl;
        exit(1);
    }

    mmf.close();
    return ret;
}

void SPM::Print(std::ostream &out)
{
    SPM::PntIter p, p_start, p_end;

    p_start = PointsBegin();
    p_end = PointsEnd();
    for (p = p_start; p != p_end; ++p)
        out << " " << (*p);
        
    out << std::endl;
}

void SPM::PrintRows(std::ostream &out) {
    std::cout << "Row Ptr: ";
    for (long unsigned i = 0; i < rowptr_size_; i++)
        std::cout << rowptr_[i] << " ";
    std::cout << std::endl;
}

void SPM::PrintElems(std::ostream &out)
{
    SPM::PntIter p, p_start, p_end;
    static int cnt = 1;

    p_start = PointsBegin();
    p_end = PointsEnd();
    for (p = p_start; p != p_end; ++p) {
        if ((*p).pattern == NULL) {
            out << std::setiosflags(std::ios::scientific)
                << row_start_ + (*p).y << " "
                << (*p).x << " " << (*p).val << " cnt:" << cnt++ << "\n";
            continue;
        }

        boost::function<void (CooElem &p)> xform_fn, rxform_fn;

        DeltaRLE *pat;
        CooElem start;
        double *vals;

        pat = (*p).pattern;
        start = static_cast<CooElem>(*p);
        vals = start.vals;
        //out << SpmTypesNames[pat->GetType()] << std::endl;
        xform_fn = GetTransformFn(type_, pat->GetType());
        rxform_fn = GetTransformFn(pat->GetType(), type_);
        if (xform_fn)
            xform_fn(start);

        DeltaRLE::Generator *g = (*p).pattern->generator(start);

        while (!g->IsEmpty()) {
            CooElem e = g->Next();
            if (rxform_fn)
                rxform_fn(e);
                
            out << std::setiosflags(std::ios::scientific) << row_start_ + e.y
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

    nr_rows_with_patterns = GetNrRows();
    nr_patterns = 0;
    nr_nzeros_block = 0;
    nr_transitions = 0;
    memset(nr_xform_patterns, 0, sizeof(nr_xform_patterns));
    for (uint64_t i = 0; i < GetNrRows(); i++) {
        uint64_t pt_size, pt_size_before;
        SpmIterOrder pt_type, pt_type_before;

        nr_patterns_before = nr_patterns;

        const SpmRowElem *elem = RowBegin(i);

        if (elem->pattern) {
            pt_size_before = pt_size = elem->pattern->GetSize();
            pt_type_before = pt_type = elem->pattern->GetType();
        } else {
            pt_size_before = pt_size = 0;
            pt_type_before = pt_type = NONE;
        }

        for (; elem != RowEnd(i); elem++) {
            if (elem->pattern) {
                ++nr_patterns;
                pt_size = elem->pattern->GetSize();
                pt_type = elem->pattern->GetType();
                nr_nzeros_block += pt_size;
                if (pt_type != pt_type_before ||
                    (pt_size_before && pt_size != pt_size_before))
                    ++nr_transitions;
                    
                ++nr_xform_patterns[elem->pattern->GetType()];
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
        ridx = GetNrRows();
        
    return PntIter(this, ridx);
}

static inline bool elem_cmp_less(const SpmCooElem &e0, const SpmCooElem &e1)
{
    int ret;

    ret = CooCmp(static_cast<CooElem>(e0), static_cast<CooElem>(e1));
    return (ret < 0);
}

static inline bool is_row_block(SpmIterOrder t)
{
    if (t > BLOCK_TYPE_START && t < BLOCK_COL_START)
        return true;
    return false;
}

static inline bool is_col_block(SpmIterOrder t)
{
    if (t > BLOCK_COL_START && t < BLOCK_TYPE_END)
        return true;
    return false;
}

static inline int gcd(int i, int j)
{
    if (j == 0)
        return i;
    return gcd(j, i%j);
}

static inline int lcm(int i, int j)
{
    if (i >= j)
        return i * j / gcd(i,j);
    else
        return i * j / gcd(j,i);
}

uint64_t SPM::FindNewRowptrSize(SpmIterOrder t)
{
    if (t == HORIZONTAL)
        return nr_rows_ + 1;
    else if (t == VERTICAL)
        return nr_cols_ + 1;
    else if (t == DIAGONAL || REV_DIAGONAL)
        return nr_rows_ + nr_cols_ + 1;
    else if (is_row_block(t))
        return nr_rows_ / IsBlockType(t) + 2;
    else if (is_col_block(t))
        return nr_cols_ / IsBlockType(t) + 2;
    else
        return 0;
}

void SPM::Transform(SpmIterOrder t, uint64_t rs, uint64_t re)
{
    PntIter p0, pe, p;
    std::vector<SpmCooElem> elems;
    std::vector<SpmCooElem>::iterator e0, ee, es;
    boost::function<void (CooElem &p)> xform_fn;

    if (type_ == t)
        return;

    xform_fn = GetTransformFn(type_, t);
    elems.reserve(elems_size_);
    
    p0 = PointsBegin(rs);
    pe = PointsEnd(re);
    for (p = p0; p != pe; ++p) {
        SpmCooElem p_new = SpmCooElem(*p);
        xform_fn(p_new);
        elems.push_back(p_new);
    }
    
    if (((type_ == HORIZONTAL || is_row_block(type_)) &&
         (t == HORIZONTAL || is_row_block(t))) ||
        ((type_ == VERTICAL || is_col_block(type_)) &&
         (t == VERTICAL || is_col_block(t)))) {
        int old_block_align, new_block_align, k;

        if (type_ == HORIZONTAL || type_ == VERTICAL)
            old_block_align = 1;
        else
            old_block_align = IsBlockType(type_);

        if (t == HORIZONTAL || t == VERTICAL)
            new_block_align = 1;
        else
            new_block_align = IsBlockType(t);

        k = lcm(old_block_align, new_block_align);
        k /= old_block_align;
        re = GetNrRows();
        
        e0 = elems.begin();
        ee = e0;
        for (uint64_t i = k; i < re; i += k) {
            es = ee;
            ee += rowptr_[i] - rowptr_[i-k];
            sort(es, ee, elem_cmp_less);
        }
        es = ee;
        ee = elems.end();
        sort(es, ee, elem_cmp_less);
    } else {
        e0 = elems.begin();
        ee = elems.end();
        sort(e0, ee, elem_cmp_less);
    }

    SetElems(e0, ee, rs + 1, 0, elems_size_, FindNewRowptrSize(t));
    elems.clear();
    type_ = t;
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
        ret = bll::bind(pnt_rmap_D, bll::_1, bll::_1, nr_rows_);
        break;
    case REV_DIAGONAL:
        ret = bll::bind(pnt_rmap_rD, bll::_1, bll::_1, nr_cols_);
        break;
    case BLOCK_ROW_TYPE_NAME(1):
        ret = bll::bind(BLOCK_ROW_RMAP_NAME(1), bll::_1, bll::_1);
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
    case BLOCK_COL_TYPE_NAME(1):
        ret = bll::bind(BLOCK_COL_RMAP_NAME(1), bll::_1, bll::_1);
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
        ret = bll::bind(pnt_map_D, bll::_1, bll::_1, nr_rows_);
        break;
    case REV_DIAGONAL:
        ret = bll::bind(pnt_map_rD, bll::_1, bll::_1, nr_cols_);
        break;
    case BLOCK_ROW_TYPE_NAME(1):
        ret = bll::bind(BLOCK_ROW_MAP_NAME(1), bll::_1, bll::_1);
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
    case BLOCK_COL_TYPE_NAME(1):
        ret = bll::bind(BLOCK_COL_MAP_NAME(1), bll::_1, bll::_1);
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
    if (rs + length > rowptr_size_ - 1)
        length = rowptr_size_ - rs - 1;

    SPM *ret = new SPM();
    uint64_t es = rowptr_[rs];
    uint64_t ee = rowptr_[rs+length];

    ret->rowptr_ = (uint64_t *) malloc((length+1) * sizeof(uint64_t));
    ret->rowptr_size_ = length + 1;
    for (uint64_t i = 0; i < ret->rowptr_size_; i++)
        ret->rowptr_[i] = rowptr_[rs+i] - es;

    ret->elems_ = &elems_[es];
    ret->elems_size_ = ee - es;
    ret->nr_rows_ = ret->rowptr_size_ - 1;
    ret->nr_cols_ = nr_cols_;
    ret->nr_nzeros_ = ret->elems_size_;
    ret->row_start_ = row_start_ + rs;
    ret->type_ = type_;
    ret->elems_mapped_ = true;
    assert(ret->rowptr_[ret->rowptr_size_-1] == ret->elems_size_);
    return ret;
}

SPM *SPM::ExtractWindow(uint64_t rs, uint64_t length)
{
    std::vector<SpmCooElem> elems;
    std::vector<SpmCooElem>::iterator elem_begin, elem_end;
    PntIter p0, pe, p;
    SPM *ret = new SPM();

    if (rs + length > rowptr_size_ - 1)
        length = rowptr_size_ - rs - 1;

    uint64_t es = rowptr_[rs];
    uint64_t ee = rowptr_[rs+length];

    assert(es <= ee);
    elems.reserve(ee - es);
    p0 = PointsBegin(rs);
    pe = PointsEnd(rs + length - 1);
    for (p = p0; p != pe; ++p)
        elems.push_back(*p);

    elem_begin = elems.begin();
    elem_end = elems.end();
    ret->SetElems(elem_begin, elem_end, rs + 1, elems_size_);
    elems.clear();
    ret->nr_rows_ = ret->rowptr_size_ - 1;
    ret->nr_cols_ = nr_cols_;
    ret->nr_nzeros_ = ret->elems_size_;
    ret->row_start_ = row_start_ + rs;
    ret->type_ = type_;
    ret->elems_mapped_ = true;
    assert(ret->rowptr_[ret->rowptr_size_-1] == ret->elems_size_);
    return ret;
}

void SPM::PutWindow(const SPM *window)
{
    assert(window);
    assert(type_ == window->type_);

    uint64_t rs = window->row_start_ - row_start_;
    uint64_t es = rowptr_[rs];

    for (uint64_t i = 0; i < window->rowptr_size_; i++)
        rowptr_[rs+i] = es + window->rowptr_[i];

    if (!window->elems_mapped_)
        memcpy(&elems_[es], window->elems_,
               window->elems_size_*sizeof(SpmRowElem));
}

SPM::Builder::Builder(SPM *spm, uint64_t nr_elems, uint64_t nr_rows) : spm_(spm)
{
    uint64_t *rowptr;
    
    if (spm_->elems_mapped_) {
        da_elems_ = dynarray_init_frombuff(sizeof(SpmRowElem),
                                           spm_->elems_size_,
                                           spm_->elems_,
                                           spm_->elems_size_);
        dynarray_seek(da_elems_, 0);
    } else {
        da_elems_ =
            dynarray_create(sizeof(SpmRowElem), 512, nr_elems);
    }
    da_rowptr_ = dynarray_create(sizeof(uint64_t), 512, nr_rows);
    rowptr = (uint64_t *) dynarray_alloc(da_rowptr_);
    *rowptr = 0;
}

SPM::Builder::~Builder()
{
    assert(da_elems_ == NULL && "da_elems_ not destroyed");
    assert(da_rowptr_ == NULL && "da_rowptr_ not destroyed");
}

SpmRowElem *SPM::Builder::AllocElem()
{
    if (spm_->elems_mapped_)
        assert(dynarray_size(da_elems_) < spm_->elems_size_ &&
               "out of bounds");
    return (SpmRowElem *) dynarray_alloc(da_elems_);
}

SpmRowElem *SPM::Builder::AllocElems(uint64_t nr)
{
    if (spm_->elems_mapped_)
        assert(dynarray_size(da_elems_) + nr < spm_->elems_size_ &&
               "out of bounds");
        
    return (SpmRowElem *) dynarray_alloc_nr(da_elems_, nr);
}

uint64_t SPM::Builder::GetElemsCnt()
{
    return dynarray_size(da_elems_);
}

void SPM::Builder::NewRow(uint64_t rdiff)
{
    uint64_t *rowptr;
    uint64_t elems_cnt;

    elems_cnt = GetElemsCnt();
    rowptr = (uint64_t *) dynarray_alloc_nr(da_rowptr_, rdiff);
    for (uint64_t i = 0; i < rdiff; i++)
        rowptr[i] = elems_cnt;
}

void SPM::Builder::Finalize()
{
    uint64_t *last_rowptr;

    last_rowptr = (uint64_t *) dynarray_get_last(da_rowptr_);
    if (*last_rowptr != dynarray_size(da_elems_))
        NewRow();

    if (!spm_->elems_mapped_ && spm_->elems_)
        free(spm_->elems_);

    if (spm_->rowptr_)
        free(spm_->rowptr_);

    if (spm_->elems_mapped_)
        assert(spm_->elems_size_ == dynarray_size(da_elems_));

    spm_->elems_size_ = dynarray_size(da_elems_);
    
    if (spm_->elems_mapped_)
        free(da_elems_);
    else
        spm_->elems_ = (SpmRowElem *) dynarray_destroy(da_elems_);

    spm_->rowptr_size_ = dynarray_size(da_rowptr_);
    spm_->rowptr_ = (uint64_t *) dynarray_destroy(da_rowptr_);
    da_elems_ = da_rowptr_ = NULL;
}

SPM::PntIter::PntIter(SPM *spm, uint64_t r_idx) : spm_(spm), row_idx_(r_idx)
{
    uint64_t *rp = spm_->rowptr_;
    uint64_t rp_size = spm_->rowptr_size_;
    
    assert(r_idx < rp_size);
    while (r_idx + 1 < rp_size && rp[r_idx] == rp[r_idx+1])
        r_idx++;
        
    row_idx_ = r_idx;
    elem_idx_ = rp[r_idx];
}

bool SPM::PntIter::operator==(const PntIter &pi)
{
    return (spm_ = pi.spm_) && (row_idx_ == pi.row_idx_) &&
           (elem_idx_ == pi.elem_idx_);
}

bool SPM::PntIter::operator!=(const PntIter &pi)
{
    return !(*this == pi);
}

void SPM::PntIter::operator++()
{
    uint64_t *rp = spm_->rowptr_;
    uint64_t rp_size = spm_->rowptr_size_;
    
    assert(elem_idx_ < spm_->elems_size_);
    assert(row_idx_ < rp_size);
    elem_idx_++;
    while (row_idx_ + 1 < rp_size && rp[row_idx_+1] == elem_idx_)
        row_idx_++;
}

SpmCooElem SPM::PntIter::operator*()
{
    SpmCooElem ret;
    SpmRowElem *e;
    DeltaRLE *p;
    
    ret.y = row_idx_ + 1;
    e = spm_->elems_ + elem_idx_;
    ret.x = e->x;
    ret.val = e->val;
    p = e->pattern;
    ret.pattern = (p == NULL) ? NULL : p->Clone();
    ret.in_pattern = e->in_pattern;
    ret.pattern_start = e->pattern_start;
    if (p != NULL)
        delete p;
        
    return ret;
}

std::ostream &SPM::PntIter::operator<<(std::ostream &out)
{
    out << "<" << std::setw(2) << row_idx_ << "," << std::setw(2)
        << elem_idx_ << ">";
    return out;
}

SPMSym *SPMSym::LoadMMF_mt(const char *mmf_file, const long nr)
{
    SPMSym *ret;
    std::ifstream mmf;

    mmf.open(mmf_file);

    if (mmf.good()) {
        ret = LoadMMF_mt(mmf, nr);
    } else {
        std::cerr << "File error!" << std::endl;
        exit(1);
    }

    mmf.close();
    return ret;
}

SPMSym *SPMSym::LoadMMF_mt(std::istream &in, const long nr)
{
    MMF mmf(in);
    return LoadMMF_mt(mmf, nr);
}

SPMSym *SPMSym::LoadMMF_mt(MMF &mmf, const long nr)
{
    SPMSym *ret, *spm_sym;
    long limit, cnt, row_start, nr_nzeros, n, nnz;
    MMF::iterator iter = mmf.begin();
    MMF::iterator iter_end = mmf.end();

    assert(mmf.GetNrRows() == mmf.GetNrCols());
    nr_nzeros = (mmf.GetNrNonzeros() + mmf.GetNrCols()) / 2;
    n = mmf.GetNrCols();
    ret = new SPMSym[nr];
    row_start = limit = cnt = 0;
    for (long i = 0; i < nr; ++i) {
        spm_sym = ret + i;
        limit = (nr_nzeros - cnt) / (nr - i);
        nnz = spm_sym->SetElems(iter, iter_end, row_start + 1, limit,
                                limit + 2 * mmf.GetNrRows() - 1,
                                mmf.GetNrRows() + 1);
        spm_sym->lower_matrix_->SetNrNonzeros(nnz - spm_sym->GetDiagonalSize());
        spm_sym->lower_matrix_->SetNrRows(spm_sym->lower_matrix_->GetNrRows());
        spm_sym->lower_matrix_->SetNrCols(n);
        spm_sym->lower_matrix_->SetRowStart(row_start);
        row_start += spm_sym->GetDiagonalSize();
        spm_sym->lower_matrix_->SetType(HORIZONTAL);
        cnt += nnz;
    }

    assert((uint64_t) cnt == (uint64_t) nr_nzeros);
    return ret;
}

void SPMSym::DivideMatrix()
{
    SPM *matrix = lower_matrix_;
    uint64_t *rowptr = matrix->GetRowPtr();
    SpmRowElem *elems = matrix->GetElems();
    uint64_t nr_nzeros = matrix->GetNrNonzeros();
    uint64_t row_start = matrix->GetRowStart();
    uint64_t nr_rows = matrix->GetNrRows();
    uint64_t nr_cols = matrix->GetNrCols();
    uint64_t rows1 = 0;
    uint64_t rows2 = 0;
    
    SPM::Builder *spmbld1 = new SPM::Builder(m1_, nr_nzeros, nr_rows + 1);
    SPM::Builder *spmbld2 = new SPM::Builder(m2_, nr_nzeros, nr_rows + 1);
    
    SpmRowElem *elem;
    
    m1_->SetType(HORIZONTAL);
    m1_->SetRowStart(row_start);
    m1_->SetNrCols(nr_cols);
    m1_->SetNrNonzeros(0);
    
    m2_->SetType(HORIZONTAL);
    m2_->SetRowStart(row_start);
    m2_->SetNrCols(nr_cols);
    m2_->SetNrNonzeros(0);
    
    for (unsigned long i = 0; i < nr_rows; i++) {
        for (unsigned long j = rowptr[i]; j < rowptr[i+1]; j++) {
            elem = elems + j;
            if (elem->x < row_start + 1) {
                if (rows1 < i) {
                    spmbld1->NewRow(i - rows1);
                    rows1 = i;
                }
                m1_->SetNrNonzeros(m1_->GetNrNonzeros() + 1);
                MakeRowElem(*elem, spmbld1->AllocElem());
            } else {
                if (rows2 < i) {
                    spmbld2->NewRow(i - rows2);
                    rows2 = i;
                }
                m2_->SetNrNonzeros(m2_->GetNrNonzeros() + 1);
                MakeRowElem(*elem, spmbld2->AllocElem());
            }
        }
    }
    
    spmbld1->Finalize();
    spmbld2->Finalize();
    
    m1_->SetNrRows(m1_->GetNrRows());
    m2_->SetNrRows(m2_->GetNrRows());
    
    delete spmbld1;
    delete spmbld2;
}

void SPMSym::MergeMatrix()
{
    SPM *m1 = GetFirstMatrix();
    SPM *m2 = GetSecondMatrix();
    SPM *temp = new SPM();
    
    SPM *matrix = lower_matrix_;
    uint64_t *rowptr1 = m1->GetRowPtr();
    uint64_t *rowptr2 = m2->GetRowPtr();
    SpmRowElem *elems1 = m1->GetElems();
    SpmRowElem *elems2 = m2->GetElems();
    uint64_t row_start = matrix->GetRowStart();
    uint64_t nr_rows = matrix->GetNrRows();
    uint64_t nr_cols = matrix->GetNrCols();
    uint64_t nr_nzeros = matrix->GetNrNonzeros();
    
    SPM::Builder *spmbld = new SPM::Builder(temp, nr_nzeros, nr_rows + 1);
    
    SpmRowElem *elem;
    
    temp->SetType(HORIZONTAL);
    temp->SetRowStart(row_start);
    temp->SetNrCols(nr_cols);
    temp->SetNrNonzeros(matrix->GetNrNonzeros());
    
    for (unsigned long i = 0; i < nr_rows; i++) {
        if (m1_->GetNrRows() > i) {
            for (unsigned long j = rowptr1[i]; j < rowptr1[i+1]; j++) {
                elem = elems1 + j;
                MakeRowElem(*elem, spmbld->AllocElem());
            }
        }
        if (m2_->GetNrRows() > i) {
            for (unsigned long j = rowptr2[i]; j < rowptr2[i+1]; j++) {
                elem = elems2 + j;
                MakeRowElem(*elem, spmbld->AllocElem());
            }
        }
        spmbld->NewRow();
    }
    
    spmbld->Finalize();
    
    temp->SetNrRows(temp->GetNrRows());
    
    delete spmbld;
    
    delete lower_matrix_;
    lower_matrix_ = temp;
    
    delete m1_;
    delete m2_;
}

void SPMSym::PrintDiagElems(std::ostream &out)
{
    uint64_t row_start = lower_matrix_->GetRowStart();
    
    for (long unsigned i = 0; i < diagonal_size_; i++)
        std::cout << row_start + i + 1 << " " << row_start + i + 1 << " "
                  << diagonal_[i] << " cnt:" << i << std::endl;
}

SPMSym::Builder::Builder(SPMSym *spm_sym, uint64_t nr_elems, uint64_t nr_rows) 
                        : spm_sym_(spm_sym)
{
    SPM *spm = spm_sym_->GetLowerMatrix();
    
    spm_bld_ = new SPM::Builder(spm, nr_elems, nr_rows);
    da_dvalues_ = dynarray_create(sizeof(double), 512, nr_rows);
}

SPMSym::Builder::~Builder()
{
    delete spm_bld_;
    assert(da_dvalues_ == NULL && "da_dvalues_ not destroyed");
}

double *SPMSym::Builder::AllocDiagElem()
{
    // assert(dynarray_size(da_dvalues_) <
    //       spm_sym_->GetLowerMatrix()->GetNrRows() && "out of bounds");

    return (double *) dynarray_alloc(da_dvalues_);
}

uint64_t SPMSym::Builder::GetDiagElemsCnt()
{
    return dynarray_size(da_dvalues_);
}

uint64_t SPMSym::Builder::GetElemsCnt()
{
    return GetDiagElemsCnt() + spm_bld_->GetElemsCnt();
}

void SPMSym::Builder::Finalize()
{
    spm_bld_->Finalize();
    spm_sym_->SetDiagonalSize(GetDiagElemsCnt());
    spm_sym_->SetDiagonal((double *) dynarray_destroy(da_dvalues_));
    da_dvalues_ = NULL;
}

void PrintTransform(uint64_t y, uint64_t x, TransformFn xform_fn,
                    std::ostream &out)
{
    uint64_t i,j;
    CooElem p;

    for (i = 1; i <= y; i++) {
        for (j = 1; j <= x;  j++) {
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

void TestTransform(uint64_t y, uint64_t x, TransformFn xform_fn,
                   TransformFn rxform_fn)
{
    uint64_t i,j;
    CooElem p0, p1;
    
    for (i = 1; i <= y; i++) {
        for (j = 1; j <= x;  j++){
            p0.y = i;
            p0.x = j;
            xform_fn(p0);
            p1.y = p0.y;
            p1.x = p0.x;
            rxform_fn(p1);
            if (p1.y != i || p1.x != j) {
                std::cout << "Error for " << i << "," << j << std::endl;
                std::cout << "Transformed: " << p0.y << "," << p0.x
                          << std::endl;
                std::cout << "rTransformed:" << p1.y << "," << p1.x
                          << std::endl;
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

	drle = RLEncode(delta);
	FOREACH(RLE<int> &v, drle){
		std::cout << "(" << v.val << "," << v.freq << ") ";
	}
	std::cout << std::endl;
}

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
