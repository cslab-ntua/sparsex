/* -*- C++ -*-
 *
 * csx.h -- Implements the CSX manager class, which responsible for creating
 *          the final CSX matrix.
 *
 * Copyright (C) 2010, Computing Systems Laboratory (CSLab), NTUA.
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_H__
#define CSX_H__

#include "ctl_ll.h"

///< CSX matrix format
typedef struct {
    uint64_t nnz, ncols, nrows, ctl_size, row_start;
    double *values;
    uint8_t *ctl;
} csx_double_t;

#include "dynarray.h"

#define PID_DELTA_BASE 0

#ifdef __cplusplus  // This file is also included from csx_llvm_tmpl.c.

#include "spm.h"
namespace csx {

/**
 *  This class is responsible for transforming the matrix from our internal
 *  representation (SPM object) to the final CSX form.
 *
 *  @see SPM
 */
class CsxManager
{
public:
    /**
     *  Pattern information that is essential for the construction of CSX.
     */
    struct PatInfo {
        PatInfo(uint8_t flag_, uint64_t nr_): flag(flag_), nr(nr_) {}
        PatInfo(): flag(0), nr(0) {}

        uint8_t flag; 	///< A unique CSX ID assigned to this pattern.
        uint64_t nr;  	///< Number of non-zero elements of this pattern.
    };
    
    typedef std::map<long,PatInfo> PatMap;
    PatMap patterns;    ///< Patterns found in matrix.
    
    CsxManager(SPM *spm)
        : spm_(spm),
          flag_avail_(0),
          row_jmps_(false),
          ctl_da_(NULL),
          last_col_(0), empty_rows_(0) {}
    
    /**
     *  Get a unique CSX ID for the pattern with SPM ID <tt>pattern_id</tt> and
     *  updates statistics for this pattern.
     *  
     *  @param pattern_id the pattern ID in the SPM format.
     *  @param nnz        number of elements included in the specific pattern.
     *  @return           the CSX pattern ID for the specified SPM pattern.
     *  @see SPM
     */
    uint8_t GetFlag(long pattern_id, uint64_t nnz);
    
    /**
     *  Transform the matrix owned by this manager into CSX form.
     *
     *  @return a handle to the newly created CSX matrix.
     */
    csx_double_t *MakeCsx();
	
    /**
     *  Checks whether row jumps exist in the matrix to be encoded in CSX
     *  format.
     *
     *  @return <tt>true</tt> if row jumps exist.
     */
    bool HasRowJmps()
    {
        return row_jmps_;
    }
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

    SPM *spm_;
    uint8_t flag_avail_;    ///< Current available flags for pattern id mapping.
    bool row_jmps_;	        ///< Whether or not row jumps included.
    
    double *values_;
    uint64_t values_idx_;
    
    dynarray_t *ctl_da_;
    uint64_t last_col_;
    bool new_row_;          ///< Marker of new_row.
    uint64_t empty_rows_;   ///< Number of empty rows since last non-empty row.
};
}   // end of csx namespace

#endif  // __cplusplus
#endif  // CSX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
