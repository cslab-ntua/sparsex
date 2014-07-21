/*
 * Copyright (C) 2011-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2014, Vasileios Karakasis
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */ 

/**
 * \file delta_sym_tmpl.c
 * \brief Template for runtime code generation.
 * 
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

static inline spx_value_t delta${bits}_case(
    uint8_t **ctl, uint8_t size, spx_value_t **values,
    spx_value_t *x, spx_value_t *y, spx_value_t *cur,
    spx_index_t *x_indx, spx_index_t *y_indx, spx_value_t scale_f)
{
	register spx_value_t yr = 0;
	register spx_value_t values_ = **values;
	register spx_index_t x_indx_ = *x_indx;
	register spx_value_t *x_ = x + x_indx_;
	register spx_value_t *ry_ = cur + x_indx_;
	register spx_index_t y_indx_ = *y_indx;
	register spx_value_t rx_ = x[y_indx_];
    register spx_index_t jmp;

    ${align_ctl}
	yr += (*x_) * values_;
	(*ry_) += rx_ * values_ * scale_f;
	(*values)++;
	for (uint8_t i = 1; i < size; i++) {
		jmp = u${bits}_get(ctl);
		x_ += jmp;
		ry_ += jmp;
		(*x_indx) += jmp;
		values_ = **values;
		yr += (*x_) * values_;
		(*ry_) += rx_ * values_ * scale_f;
		(*values)++;
	}

	return yr * scale_f;
}
