/* -*- C++ -*-
 *
 * csx.h -- The CSX Manager interface
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_H__
#define CSX_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "ctl_ll.h"

#ifdef __cplusplus
}
#endif

///< CSX matrix format
typedef struct {
    uint64_t nnz, ncols, nrows, ctl_size, row_start;
    double *values;
    uint8_t *ctl;
} csx_double_t;

typedef struct {
    csx_double_t *lower_matrix;
    double *dvalues;
} csx_double_sym_t;

///> Deletes all the fields of CSX format.
void DestroyCsx(csx_double_t *csx);
    
///> Deletes all the fields of CSX-Sym format.
void DestroyCsxSym(csx_double_sym_t *csx_sym);

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
        PatInfo(uint8_t flag_, uint64_t npatterns_, uint64_t nr_): 
                flag(flag_), npatterns(npatterns_), nr(nr_) {}
        PatInfo(): flag(0), npatterns(0), nr(0) {}

        uint8_t flag;       // A unique CSX ID assigned to this pattern.
        uint64_t npatterns; // Number of patterns of this kind.
        uint64_t nr;        // Number of non-zero elements of this pattern.
    };

    typedef std::map<long,PatInfo> PatMap;
    PatMap patterns;        // Patterns found in matrix.

    CsxManager(SPM *spm)
        : spm_(spm),
          flag_avail_(0),
          row_jmps_(false),
          full_column_indices_(false),
          ctl_da_(NULL),
          last_col_(0), empty_rows_(0) {}
    
    CsxManager(SPMSym *spm_sym)
        : spm_sym_(spm_sym),
          flag_avail_(0),
          row_jmps_(false),
          full_column_indices_(false),
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
     *  @param symmetric determines if the matrix is symmetric (only the lower
     *                   triangle half must be encoded) or not.
     *  @return          a handle to the newly created CSX matrix or to CSX-Sym
     *                   lower triangle half part of the matrix.
     */
    csx_double_t *MakeCsx(bool symmetric);
    
    /**
     *  Transform the matrix owned by this manager into CSX form.
     *
     *  @return a handle to the newly created CSX-Sym matrix.
     */
    csx_double_sym_t *MakeCsxSym();

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

    bool HasFullColumnIndices()
    {
        return full_column_indices_;
    }

    void SetFullColumnIndices(bool val)
    {
        full_column_indices_ = val;
    }

private:
    /**
     *  Transform a row of the matrix into CSX form.
     *
     *  @rstart first element of the row.
     *  @rend   last element of the row.
     */
    //void DoRow(const SpmRowElem *rstart, const SpmRowElem *rend);
    //void DoSymRow(const SpmRowElem *rstart, const SpmRowElem *rend);
    void DoRow(const SpmElem *rstart, const SpmElem *rend);
    void DoSymRow(const SpmElem *rstart, const SpmElem *rend);

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
     */
    //void AddPattern(const SpmRowElem &elem);
    void AddPattern(const SpmElem &elem);

    /**
     *  @param xs  elements found before pattern.
     */
    //void PreparePat(std::vector<uint64_t> &xs, const SpmRowElem &elem);
    void PreparePat(std::vector<uint64_t> &xs, const SpmElem &elem);

    SPM *spm_;
    SPMSym *spm_sym_;
    uint8_t flag_avail_;     // Current available flags for pattern id mapping.
    bool row_jmps_;          // Whether or not row jumps included.
    bool full_column_indices_;  // use full 32-bit indices instead of deltas.
                                
    double *values_;
    uint64_t values_idx_;

    dynarray_t *ctl_da_;
    uint64_t last_col_;
    bool new_row_;          // Marker of new_row.
    uint64_t empty_rows_;   // Number of empty rows since last non-empty row.
};
}   // end of csx namespace

#endif  // __cplusplus
#endif  // CSX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
