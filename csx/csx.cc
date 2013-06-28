/*
 * csx.cc --  CSX utilities
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "csx.h"
#include "ctl_ll.h"

void DestroyCsx(csx_double_t *csx)
{
#ifdef SPM_NUMA
    free_interleaved(csx->ctl, csx->ctl_size*sizeof(*csx->ctl));
    free_interleaved(csx->values, csx->nnz*sizeof(*csx->values));
    free_interleaved(csx, sizeof(*csx));
#else
    free(csx->ctl);
    free(csx->values);
    free(csx);
#endif
}

void DestroyCsxSym(csx_double_sym_t *csx_sym)
{
#ifdef SPM_NUMA
    uint64_t diag_size = csx_sym->lower_matrix->nrows;
#endif    
        DestroyCsx(csx_sym->lower_matrix);
#ifdef SPM_NUMA
    numa_free(csx_sym->dvalues, diag_size*sizeof(*csx_sym->dvalues));
    numa_free(csx_sym, sizeof(*csx_sym));
#else
    free(csx_sym->dvalues);
    free(csx_sym);
#endif
}

#if 0

void PrintCsx(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;

    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        DoPrint(spm_mt->spm_threads[i].spm);
    }
}

void PrintElement(uint64_t row, uint64_t col, double val)
{
    std::cout << "(" << std::setw(2) << row << "," << std::setw(2) << col << "): "
              << val << std::endl;
}

void DoPrint(void *spm)
{
    csx_double_t *csx = (csx_double_t *) spm;
	double *v = csx->values;
	uint8_t *ctl = csx->ctl;
	uint8_t *ctl_end = ctl + csx->ctl_size;
	uint8_t size, flags, patt_id, type;
    uint8_t delta, row_dim, col_dim;
    uint64_t cur_row = csx->row_start;
    uint64_t ucol = 0;
    uint64_t i_end;

	do {
        flags = *ctl++;
		size = *ctl++;
		if (test_bit(&flags, CTL_NR_BIT)) {
            if (test_bit(&flags, CTL_RJMP_BIT)) {
                cur_row += ul_get(&ctl);
            } else {
                cur_row++;
            }
            ucol = 0;
		}
        if (false)//(full_column_indices)
            ucol = u32_get(&ctl);// + last_col;
        else
            ucol += ul_get(&ctl);

		patt_id = flags & CTL_PATTERN_MASK;
        type = csx->id_map[patt_id] / CSX_PID_OFFSET;
        delta = csx->id_map[patt_id] % CSX_PID_OFFSET;

        if (type >= 6 && type <= 13) {
            row_dim = type - BLOCK_TYPE_START;
            col_dim = delta;
        } else if (type >= 15 && type <= 22) {
            col_dim = type - BLOCK_COL_START;
            row_dim = delta;
        }

        std::cout << "New encoded element:\n";
        printf("type:%"PRIu8 " ", type);
        printf("delta:%"PRIu8 " ", delta);
        printf("size:%"PRIu8 " ", size);
        std::cout << std::endl;

        switch (type) {
        case 0:  // DELTAS
            PrintElement(cur_row + 1, ucol + 1, *v);
            v++;
            if (delta == 8) {
                for (uint8_t i = 1; i < size; i++) {
                    ucol += u8_get(&ctl);
                    PrintElement(cur_row + 1, ucol + 1, *v);
                    v++;
                }
            } else if (delta == 16) {
                align_ptr(&ctl, 2);
                for (uint8_t i = 1; i < size; i++) {
                    ucol += u16_get(&ctl);
                    PrintElement(cur_row + 1, ucol + 1, *v);
                    v++;
                }
            } else if (delta == 32){
                align_ptr(&ctl, 4);
                for (uint8_t i = 1; i < size; i++) {
                    ucol += u32_get(&ctl);
                    PrintElement(cur_row + 1, ucol + 1, *v);
                    v++;
                }
            } else {
                align_ptr(&ctl, 8);
                for (uint8_t i = 1; i < size; i++) {
                    ucol += u64_get(&ctl);
                    PrintElement(cur_row + 1, ucol + 1, *v);
                    v++;
                }
            }
            break;
        case 1: //  HORIZONTAL
            i_end = size * delta;
            for (uint64_t i = 0; i < i_end; i += delta) {
                PrintElement(cur_row + 1, ucol + i + 1, *v);
                v++;
            }
            ucol += i_end - delta;
            break;
        case 2: //VERTICAL
            i_end = size * delta;
            for (uint64_t i = 0; i < i_end; i += delta) {
                PrintElement(cur_row + i + 1, ucol + 1, *v);
                v++;
            }
            break;
        case 3: //DIAGONAL
            i_end = size * delta;
            for (uint64_t i = 0; i < i_end; i += delta) {
                PrintElement(cur_row + i + 1, ucol + i + 1, *v);
                v++;
            }
            break;
        case 4: //REV_DIAGONAL
            i_end = size * delta;
            for (uint64_t i = 0; i < i_end; i += delta) {
                PrintElement(cur_row + i + 1, ucol - i + 1, *v);
                v++;
            }
            break;
        case 6 ... 13:  //BLOCK_R1 ... BLOCK_COL_START - 1
            for (uint64_t i = 0; i < col_dim; i++) {
                for (uint64_t j = 0; j < row_dim; j++) {
                    PrintElement(cur_row + j + 1, ucol + i + 1, *v);
                    v++;
                }
            } 
            break;
        case 15 ... 22: //BLOCK_COL_START ... BLOCK_TYPE_END - 1
            for (uint64_t i = 0; i < row_dim; i++) {
                for (uint64_t j = 0; j < col_dim; j++) {
                    PrintElement(cur_row + i + 1, ucol + j + 1, *v);
                    v++;
                }
            } 
            break;
        default:
            assert(0 && "unknown pattern type");
        }
	} while (ctl < ctl_end);
}

#endif
// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
