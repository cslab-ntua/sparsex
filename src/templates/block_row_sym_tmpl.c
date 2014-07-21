/*
 * Copyright (C) 2011-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2014, Vasileios Karakasis
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */ 

/**
 * \file block_row_sym_tmpl.c
 * \brief Template for runtime code generation.
 * 
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

static inline spx_value_t block_row_${r}x${c}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t *x, spx_value_t *y, spx_value_t *cur,
    spx_index_t *x_indx, spx_index_t *y_indx, spx_value_t scale_f)
{
	register spx_value_t values_ = **values;
	register spx_index_t x_indx_ = *x_indx;
	register spx_value_t *x_ = x + x_indx_;
	register spx_value_t *ry_ = cur + x_indx_;
	register spx_index_t y_indx_ = *y_indx;
	register spx_value_t *y_ = y + y_indx_;
	register spx_value_t *rx_ = x + y_indx_;
	register spx_value_t curr_x_ = *x_;
	register spx_value_t curr_ry_ = 0;

	for (spx_index_t i = 0; i < ${c}; i++) {
		for (spx_index_t j = 0; j < ${r}; j++) {
			y_[j] += curr_x_ * values_ * scale_f;
			curr_ry_ += rx_[j] * values_;
			(*values)++;
			values_ = **values;
		}
		
		(*ry_) += curr_ry_ * scale_f;
		x_++;
		ry_++;
	    curr_x_ = *x_;
	    curr_ry_ = 0;
	}

	return 0;
}
