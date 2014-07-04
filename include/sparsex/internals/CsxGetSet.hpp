/*
 * CsxGetSet.hpp -- Get/Set routines for CSX.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_GET_SET_HPP
#define SPARSEX_INTERNALS_CSX_GET_SET_HPP

#include "sparsex/internals/Config.hpp"
#include "sparsex/internals/Csx.hpp"
#include "sparsex/internals/Encodings.hpp"
#include "sparsex/internals/Element.hpp"
#include "sparsex/internals/SpmMt.hpp"

using namespace csx;

/**
 *  Returns the value of the corresponding element (0 if element doesn't exist)
 *  Indexes row/col are assumed to be 1-based.
 */
template<typename IndexType, typename ValueType>
bool GetValueCsx(void *spm, IndexType row, IndexType col, ValueType *value);
template<typename IndexType, typename ValueType>
bool GetValueCsxSym(void *spm, IndexType row, IndexType col, ValueType *value);

/**
 *  Sets the value of the corresponding element.
 *  Returns false if the element doesn't exist.
 *  Indexes row/col are assumed to be 1-based.
 */
template<typename IndexType, typename ValueType>
bool SetValueCsx(void *spm, IndexType row, IndexType col, ValueType value);
template<typename IndexType, typename ValueType>
bool SetValueCsxSym(void *spm, IndexType row, IndexType col, ValueType value);

/* Helper-function declarations. */
template<typename IndexType, typename ValueType>
bool SearchDelta(uint8_t **ctl_ptr, ValueType **v_ptr, uint8_t delta,
                 uint8_t size, IndexType col, IndexType& cur_col,
                 ValueType& value, size_t rows_checked, bool mode);
template<typename IndexType, typename ValueType>
bool SearchHorizontal(ValueType **v_ptr, uint8_t delta, uint8_t size,
                      IndexType col, IndexType& ucol, ValueType& value,
                      size_t rows_checked, bool mode);
template<typename IndexType, typename ValueType>
bool SearchVertical(ValueType **v_ptr, uint8_t delta, uint8_t size,
                    IndexType col, IndexType row, IndexType ucol,
                    IndexType current_row, ValueType& value,
                    size_t rows_checked, bool mode);
template<typename IndexType, typename ValueType>
bool SearchDiagonal(ValueType **v_ptr, uint8_t delta, uint8_t size,
                    IndexType col, IndexType row, IndexType ucol,
                    IndexType current_row, ValueType& value,
                    size_t rows_checked, bool mode, bool rev);
template<typename IndexType, typename ValueType>
bool SearchBlock(ValueType **v_ptr, uint8_t delta, uint8_t size,
                 uint8_t type,IndexType col, IndexType row,
                 IndexType ucol, IndexType current_row, ValueType& value,
                 size_t rows_checked, bool isrow, bool mode);
template<typename IndexType, typename ValueType>
bool SearchValue(void *spm, IndexType row, IndexType col,
                 ValueType& value, bool mode);


/* Function definitions. */
template<typename IndexType, typename ValueType>
bool GetValueCsx(void *spm, IndexType row, IndexType col, ValueType *value)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    CsxMatrix<IndexType, ValueType> *csx;
    bool found = false;

    //Add boundary check for row/col
    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        csx = (CsxMatrix<IndexType, ValueType> *) spm_mt->spm_threads[i].spm;
        if ((size_t) row <= (csx->row_start + csx->nrows)) {
            found = SearchValue(csx, row, col, *value, false);
            break;
        }
    }

    return found;
}

template<typename IndexType, typename ValueType>
bool GetValueCsxSym(void *spm, IndexType row, IndexType col, ValueType *value)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    CsxSymMatrix<IndexType, ValueType> *csx_sym;
    CsxMatrix<IndexType, ValueType> *csx;
    IndexType tmp;
    bool found = false;

    //Add boundary check for row/col
    if (row < col) {
        tmp = row;
        row = col;
        col = tmp;
    }
    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        csx_sym = (CsxSymMatrix<IndexType, ValueType> *)
            spm_mt->spm_threads[i].spm;
        csx = (CsxMatrix<IndexType, ValueType> *) csx_sym->lower_matrix;
        // Check if element is on main diagonal
        if (row == col && ((size_t) row <= (csx->row_start + csx->nrows))) {
            return csx_sym->dvalues[row - csx->row_start - 1];
        } else if ((size_t) row <= (csx->row_start + csx->nrows)) {
            found = SearchValue(csx, row, col, *value, false);
            break;
        }
    }

    return found;
}

template<typename IndexType, typename ValueType>
bool SetValueCsx(void *spm, IndexType row, IndexType col, ValueType value)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    CsxMatrix<IndexType, ValueType> *csx;
    bool found = false;

    //Add boundary check for row/col
    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        csx = (CsxMatrix<IndexType, ValueType> *) spm_mt->spm_threads[i].spm;
        if ((size_t) row <= (csx->row_start + csx->nrows)) {
            found = SearchValue(csx, row, col, value, true);
            break;
        }
    }

    return found;
}

template<typename IndexType, typename ValueType>
bool SetValueCsxSym(void *spm, IndexType row, IndexType col, ValueType value)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    CsxSymMatrix<IndexType, ValueType> *csx_sym;
    CsxMatrix<IndexType, ValueType> *csx;
    IndexType tmp;
    bool found = false;

    //Add boundary check for row/col
    if (row < col) {
        tmp = row;
        row = col;
        col = tmp;
    }
    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        csx_sym =
            (CsxSymMatrix<IndexType, ValueType> *) spm_mt->spm_threads[i].spm;
        csx = (CsxMatrix<IndexType, ValueType> *) csx_sym->lower_matrix;
        // Check if element is on main diagonal
        if (row == col && ((size_t) row <= (csx->row_start + csx->nrows))) {
            csx_sym->dvalues[row - csx->row_start - 1] = value;
            return true;
        } else if ((size_t) row <= (csx->row_start + csx->nrows)) {
            found = SearchValue(csx, row, col, value, true);
            break;
        }
    }

    return found;
}

#ifdef PTR_ALIGN
#define ALIGN(buf,a) (uint8_t *) (((unsigned long) (buf) + (a-1)) & ~(a-1))
void align_ptr(uint8_t **ctl, int align)
{
    *ctl = ALIGN(*ctl, align);
}
#endif

template<typename IndexType, typename ValueType>
bool SearchValue(void *spm, IndexType row, IndexType col, ValueType& value,
                 bool mode)
{
    CsxMatrix<IndexType, ValueType> *csx =
        (CsxMatrix<IndexType, ValueType> *) spm;
    IndexType current_row = row - 1;
    row_info_t *row_info = csx->rows_info;
    uint8_t *ctl_ptr;
    uint8_t *ctl_end = csx->ctl + csx->ctl_size;
    ValueType *v_ptr;
    uint8_t size, flags;
    uint8_t patt_id, type, delta;
    IndexType ucol = 0;
    size_t cnt = 0;
    size_t rows_checked = 0;

    do {
        /* Skip row if it is empty or it fails the span check */
        if ((size_t) current_row != csx->row_start) {
            if (((size_t) (row - 1) >
                 ((size_t) current_row + row_info[current_row -
                                                  csx->row_start].span))
                || (row_info[current_row - csx->row_start].rowptr ==
                    row_info[current_row - csx->row_start - 1].rowptr)) {
                rows_checked++;
                current_row--;
                continue;
            }
        }

        /* Update pointers in ctl and values */
        ctl_ptr = &(csx->ctl[row_info[current_row - csx->row_start].rowptr]);
        v_ptr = &(csx->values[row_info[current_row - csx->row_start].valptr]);
        ucol = 0;
        cnt = 0;

        do {
            flags = *ctl_ptr++;
            size = *ctl_ptr++;
            if (test_bit(&flags, CTL_NR_BIT)) {
                if (test_bit(&flags, CTL_RJMP_BIT))
                    ul_get(&ctl_ptr);
                /* End of line */
                if (cnt != 0) break;
            }
            cnt++;
#if SPX_USE_NUMA
            ucol = u32_get(&ctl_ptr);
#else
            ucol += ul_get(&ctl_ptr);
#endif
            patt_id = flags & CTL_PATTERN_MASK;
            type = csx->id_map[patt_id] / PatternIdOffset;
            delta = csx->id_map[patt_id] % PatternIdOffset;

            switch(type) {
            case Encoding::None:
                if (SearchDelta(&ctl_ptr, &v_ptr, delta, size, col, ucol, value,
                                rows_checked, mode))
                    return true;
                break;
            case Encoding::Horizontal:
                if (SearchHorizontal(&v_ptr, delta, size, col, ucol, value,
                                     rows_checked, mode))
                    return true;
                break;
            case Encoding::Vertical:
                if (SearchVertical(&v_ptr, delta, size, col, row, ucol,
                                   current_row, value, rows_checked, mode))
                    return true;
                break;
            case Encoding::Diagonal:
                if (SearchDiagonal(&v_ptr, delta, size, col, row, ucol,
                                   current_row, value, rows_checked, mode,
                                   false))
                    return true;
                break;
            case Encoding::AntiDiagonal:
                if (SearchDiagonal(&v_ptr, delta, size, col, row, ucol,
                                   current_row, value, rows_checked, mode,
                                   true))
                    return true;
                break;
            case Encoding::BlockRow1:
            case Encoding::BlockRow2:
            case Encoding::BlockRow3:
            case Encoding::BlockRow4:
            case Encoding::BlockRow5:
            case Encoding::BlockRow6:
            case Encoding::BlockRow7:
            case Encoding::BlockRow8:
                if (SearchBlock(&v_ptr, delta, size, type, col, row, ucol,
                                current_row, value, rows_checked, true, mode))
                    return true;
                break;
            case Encoding::BlockCol1:
            case Encoding::BlockCol2:
            case Encoding::BlockCol3:
            case Encoding::BlockCol4:
            case Encoding::BlockCol5:
            case Encoding::BlockCol6:
            case Encoding::BlockCol7:
            case Encoding::BlockCol8:
                if (SearchBlock(&v_ptr, delta, size, type, col, row, ucol,
                                current_row, value, rows_checked, false, mode))
                    return true;
                break;
            default:
                assert(0 && "unknown pattern type");
            }
        } while (ctl_ptr < ctl_end);

        if ((size_t) current_row == csx->row_start)
            return false;
        rows_checked++;
        current_row--;
    } while ((size_t) current_row >= csx->row_start);

    return false;
}

template<typename ValueType>
static void GetSet(ValueType& val, ValueType *cur_val, bool mode)
{
    if (mode)
        *cur_val = val;
    else
        val = *cur_val;
}

template<typename IndexType, typename ValueType>
bool SearchDelta(uint8_t **ctl_ptr, ValueType **v_ptr, uint8_t delta,
                 uint8_t size, IndexType col, IndexType& cur_col,
                 ValueType& value, size_t rows_checked, bool mode)
{
    if (delta == 8) {
        if ((rows_checked == 0) && ((col - 1) == cur_col)) {
            GetSet(value, *v_ptr, mode);
            return true;
        }
        (*v_ptr)++;
        for (uint8_t i = 1; i < size; i++) {
            cur_col += u8_get(ctl_ptr);
            if ((rows_checked == 0) && ((col - 1) == cur_col)) {
                GetSet(value, *v_ptr, mode);
                return true;
            }
            (*v_ptr)++;
        }
    } else if (delta == 16) {
#ifdef PTR_ALIGN
        align_ptr(ctl_ptr, 2);
#endif
        if ((rows_checked == 0) && ((col - 1) == cur_col)) {
            GetSet(value, *v_ptr, mode);
            return true;
        }
        (*v_ptr)++;
        for (uint8_t i = 1; i < size; i++) {
            cur_col += u16_get(ctl_ptr);
            if ((rows_checked == 0) && ((col - 1) == cur_col)) {
                GetSet(value, *v_ptr, mode);
                return true;
            }
            (*v_ptr)++;
        }
    } else if (delta == 32){
#ifdef PTR_ALIGN
        align_ptr(ctl_ptr, 4);
#endif
        if ((rows_checked == 0) && ((col - 1) == cur_col)) {
            GetSet(value, *v_ptr, mode);
            return true;
        }
        (*v_ptr)++;
        for (uint8_t i = 1; i < size; i++) {
            cur_col += u32_get(ctl_ptr);
            if ((rows_checked == 0) && ((col - 1) == cur_col)) {
                GetSet(value, *v_ptr, mode);
                return true;
            }
            (*v_ptr)++;
        }
    } else {
#ifdef PTR_ALIGN
        align_ptr(ctl_ptr, 8);
#endif
        if ((rows_checked == 0) && ((col - 1) == cur_col)) {
            GetSet(value, *v_ptr, mode);
            return true;
        }
        (*v_ptr)++;
        for (uint8_t i = 1; i < size; i++) {
            cur_col += u64_get(ctl_ptr);
            if ((rows_checked == 0) && ((col - 1) == cur_col)) {
                GetSet(value, *v_ptr, mode);
                return true;
            }
            (*v_ptr)++;
        }
    }

    return false;
}

template<typename IndexType, typename ValueType>
bool SearchHorizontal(ValueType **v_ptr, uint8_t delta, uint8_t size,
                      IndexType col, IndexType& ucol, ValueType& value,
                      size_t rows_checked, bool mode)
{
    if ((rows_checked == 0) && ((col - 1) >= ucol && 
                                (col - 1) < (ucol + size * delta))) {
        for (IndexType i = 0; i < size * delta; i += delta) {
            if ((col - 1) == (ucol + i)) {
                GetSet(value, *v_ptr, mode);
                return true;
            }
            (*v_ptr)++;
        }
    } else {
        *v_ptr += size;
    }
    ucol += (size - 1) * delta;

    return false;
}

template<typename IndexType, typename ValueType>
bool SearchVertical(ValueType **v_ptr, uint8_t delta, uint8_t size,
                    IndexType col, IndexType row, IndexType ucol,
                    IndexType current_row, ValueType& value,
                    size_t rows_checked, bool mode)
{
    if (rows_checked == 0) {
        if ((col - 1) == ucol) {
            GetSet(value, *v_ptr, mode);
            return true;
        }
        *v_ptr += size;
    } else {
        if ((col - 1) == ucol) {
            (*v_ptr)++;
            /* Find row in pattern */
            for (IndexType i = current_row + delta;
                 i < current_row + size * delta; i += delta) {
                if (i == (row - 1)) {
                    GetSet(value, *v_ptr, mode);
                    return true;
                }
                (*v_ptr)++;
            }
        } else {
            *v_ptr += size;
        }
    }

    return false;
}

template<typename IndexType, typename ValueType>
bool SearchDiagonal(ValueType **v_ptr, uint8_t delta, uint8_t size,
                    IndexType col, IndexType row, IndexType ucol,
                    IndexType current_row, ValueType& value,
                    size_t rows_checked, bool mode, bool rev)
{
    if (rows_checked == 0) {
        if ((col - 1) == ucol) {
            GetSet(value, *v_ptr, mode);
            return true;
        }
        *v_ptr += size;
    } else {
        if (((!rev && (col - 1) >= ucol &&
              (col - 1) <= (ucol + delta * (size - 1))) ||
             (rev && (col - 1) <= ucol &&
              (col - 1) >= (ucol - delta * (size - 1)))) &&
            (row - 1) <= (current_row + delta * (size - 1)) &&
            (row - 1) > current_row) {
            (*v_ptr)++;
            ucol += delta * (2 * !rev - 1);
            /* Find row in pattern */
            for (IndexType i = current_row + delta;
                 i < current_row + size * delta; i += delta) {
                if (((i == (row - 1)) && ((col - 1) == ucol))) {
                    GetSet(value, *v_ptr, mode);
                    return true;
                }
                (*v_ptr)++;
                ucol += delta * (2 * !rev - 1);
            }
        } else {
            *v_ptr += size;
        }
    }

    return false;
}

template<typename IndexType, typename ValueType>
bool SearchBlock(ValueType **v_ptr, uint8_t delta, uint8_t size,
                 uint8_t type,IndexType col, IndexType row,
                 IndexType ucol, IndexType current_row, ValueType& value,
                 size_t rows_checked, bool isrow, bool mode)
{
    IndexType row_dim = 0, col_dim = 0;
    if (isrow) {
        row_dim = type - Encoding::BlockRowMin + 1;//CHECK BLOCK_TYPE_START
        col_dim = delta;
    } else {
        row_dim = delta;
        col_dim = type - Encoding::BlockColMin + 1;//CHAECK BLOCK_COL_START
    }

    if (rows_checked == 0) {
        if ((col - 1) >= ucol && (col - 1) < (ucol + col_dim)) {
            /* Element is somewhere in here */
            if (isrow) 
                *v_ptr += row_dim * (col - ucol - 1);
            else 
                *v_ptr += col - ucol - 1;
            GetSet(value, *v_ptr, mode);
            return true;
        } else {  
            *v_ptr += size;
        }
    } else {
        if ((col - 1) >= ucol && (col - 1) < (ucol + col_dim) &&
            (row - 1) < (current_row + row_dim) && (row - 1) > current_row) {
            /* Element is somewhere in here */
            if (isrow) {
                *v_ptr += (row - 1) - current_row + row_dim * (col - ucol - 1);
            } else { 
                *v_ptr += col_dim * ((row - 1) - current_row) + (col - ucol - 1);
            }
            GetSet(value, *v_ptr, mode);
            return true;
        } else {                    
            *v_ptr += size;
        }
    }

    return false;
}

#endif  // SPARSEX_INTERNALS_CSX_GET_SET_HPP
