/* -*- C++ -*-
 *
 * csx.h -- CSX Manager
 *
 * Copyright (C) 2010, Computing Systems Laboratory (CSLab), NTUA.
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_H__
#define CSX_H__

#include "ctl_ll.h"

typedef struct {
    uint64_t nnz, ncols, nrows, ctl_size, row_start;
    double *values;
    uint8_t *ctl;
} csx_double_t;

#ifdef __cplusplus

#include "spm.h"

extern "C" {
    #include "dynarray.h"
}

#define PID_DELTA_BASE 0

namespace csx {

/**
 *  Responsible for transforming the matrix from Spm object to its final CSX 
 *  form.
 */
class CsxManager
{
public:
    SPM *spm;
 
    /**
     *  Keeps pattern info.
     */
    class PatInfo {
    public:
        uint8_t flag; 	///< Flag allocated for this pattern.
        uint64_t nr;  	///< Number of non-zero elemenets for this pattern.
        
        PatInfo(uint8_t flag_, uint64_t nr_): flag(flag_), nr(nr_) {}
        PatInfo(): flag(0), nr(0) {}
    };
    
    typedef std::map<long,PatInfo> PatMap;
    
    PatMap patterns;    ///< Patterns found in matrix.
    
    uint8_t flag_avail; ///< Current available flags for pattern id mapping.
    bool row_jmps;	///< Whether or not row jumps included.
    
    double *values;
    uint64_t values_idx;
    
    dynarray_t *ctl_da;
    uint64_t last_col;
    bool new_row;           ///< Marker of new_row.
    uint64_t empty_rows;    ///< Number of empty rows since last non-empty row.
    
    CsxManager(SPM *spm_) : spm(spm_), flag_avail(0), row_jmps(false), 
                            ctl_da(NULL), last_col(0), empty_rows(0) {}
    
    /**
     *  Finds pattern info, update its values and return its flag.
     *  
     *  @param pattern_id unique number for every type of pattern.
     *  @param nnz        number of elements included in the specific pattern.
     *  @return           flag used for the specific type of pattern.
     */
    uint8_t GetFlag(long pattern_id, uint64_t nnz);
    
    /**
     *  Transform the matrix into CSX form.
     *
     *  @return characteristics of matrix in CSX form.
     */
    csx_double_t *MakeCsx();
	
private:
    /**
     *  Transform a row of the matrix into CSX form.
     *
     *  @rstart first element of the row.
     *  @rend   last element of the row.
     */
    void DoRow(const SpmRowElem *rstart, const SpmRowElem *rend);
    
    /**
     *  Set flags that concern change of row.
     *
     *  @flags flags that are going to be updated.
     */
    void UpdateNewRow(uint8_t *flags);
    
    /**
     *  Set all the flags of pattern info for elements that do not adhere to
     *  any type of pattern.
     *
     *  @param xs vector with the columns of elements.
     */
    void AddXs(std::vector<uint64_t> &xs);
    
    /**
     *  Set all the flags of pattern info for elements that adhere to a type of
     *  pattern.
     *
     *  @param elem elements of current pattern.
     *  @param jmp  ????????????
     */
    void AddPattern(const SpmRowElem &elem, uint64_t jmp);
    
    /**
     *  ???????????????????
     */
    uint64_t PreparePat(std::vector<uint64_t> &xs, const SpmRowElem &elem);
};

}

#endif

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
