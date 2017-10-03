/*
 * Copyright (C) 2012-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2014  Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file SparseMatrix.hpp
 * \brief Generic representation of a sparse matrix
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_SPARSE_MATRIX_HPP
#define SPARSEX_INTERNALS_SPARSE_MATRIX_HPP

#include <sparsex/internals/Csr.hpp>
#include <sparsex/internals/CsxBuild.hpp>
#include <sparsex/internals/CsxGetSet.hpp>
#include <sparsex/internals/CsxSaveRestore.hpp>
#include <sparsex/internals/CsxUtil.hpp>
#include <sparsex/internals/Mmf.hpp>
#include <sparsex/internals/Rcm.hpp>
#include <sparsex/internals/SparseInternal.hpp>
#include <sparsex/internals/TimerCollection.hpp>
#include <sparsex/internals/Types.hpp>
#include <boost/interprocess/detail/move.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/static_assert.hpp>
#include <sstream>

using namespace std;
using namespace sparsex::utilities;

namespace sparsex {
  namespace csx {

    namespace internal {

      template<class MatrixType>
      struct matrix_traits {
	typedef typename MatrixType::idx_t idx_t;
	typedef typename MatrixType::val_t val_t;
      };

      /**
       *  Helper struct that generates a compile time error if SparseMatrix
       *  class is instantiated with something other than the allowed types,
       *  i.e., an integral and a floating point type.
       */
      template<typename T1, typename T2>
      void valid_instantiation
      (typename boost::enable_if<boost::is_integral<T1> >::type *dum1 = 0,
       typename boost::enable_if<boost::is_floating_point<T2> >::type *dum2 = 0)
      {}

      template<typename index_type, typename value_type>
      struct allow_instantiation
      {
	allow_instantiation()
	{
	  valid_instantiation<index_type, value_type>();
	}
      };

      template<int T>
      struct Sym
      {
	enum { value = T };
      };

    } // end namespace internal

    template<class InputPolicy>
    class SparseMatrix : public InputPolicy
    {
    public:
      typedef typename internal::matrix_traits<InputPolicy>::idx_t idx_t;
      typedef typename internal::matrix_traits<InputPolicy>::val_t val_t;

      // CSR-specific constructor
      SparseMatrix(const idx_t *rowptr, const idx_t *colind,
		   const val_t *values,
		   idx_t nr_rows, idx_t nr_cols, bool zero_based)
        : InputPolicy(boost::interprocess::forward<idx_t*>
		      (const_cast<idx_t*>(rowptr)),
                      boost::interprocess::forward<idx_t*>
		      (const_cast<idx_t*>(colind)), 
                      boost::interprocess::forward<val_t*>
		      (const_cast<val_t*>(values)),
                      boost::interprocess::forward<idx_t>(nr_rows),
                      boost::interprocess::forward<idx_t>(nr_cols),
                      boost::interprocess::forward<bool>(zero_based)),
          csx_(0)
      {
        DoCreateTimers();
      }

      // MMF-specific constructor
      SparseMatrix(const char *filename)
        : InputPolicy(boost::interprocess::forward<const char*>(filename)),
          csx_(0)
      {
        DoCreateTimers();
      }

      ~SparseMatrix()
      {
        // if (csx_)
        //     Destroy();
      }

      void Reorder(vector<size_t>& perm)
      {
        DoReorder_RCM(*this, perm);
      }

      spm_mt_t *CreateCsx()
      {
        spm_mt_t *spm = 0;
        RtConfig &config = RtConfig::GetInstance();

        try {
	  if (config.GetProperty<bool>(RtConfig::MatrixSymmetric)) {
	    LOG_INFO << "Format: CSX-sym\n";
	    spm = DoCreateCsx(internal::Sym<true>());
	  } else {
	    LOG_INFO << "Format: CSX\n";
	    spm = DoCreateCsx(internal::Sym<false>());
	  }
        } catch (std::exception &e) {
	  LOG_ERROR << e.what() << "\n";
	  exit(1);
        }

        return spm;
      }

      val_t GetValue(idx_t row_idx, idx_t col_idx)
      {
        val_t value = 0;

        if (!csx_) {
	  LOG_ERROR << "matrix hasn't been encoded yet\n";
	  exit(1);
        } else {
	  if (csx_->symmetric) {
	    if (!GetValueCsxSym<idx_t, val_t>(csx_, row_idx,
					      col_idx, &value))
	      LOG_WARNING << "matrix entry not found\n";
	  } else {
	    if (!GetValueCsx<idx_t, val_t>(csx_, row_idx,
					   col_idx, &value))
	      LOG_WARNING << "matrix entry not found\n";
	  }
        }

        return value;
      }

      bool SetValue(idx_t row_idx, idx_t col_idx, val_t value)
      {
        if (!csx_) {
	  LOG_ERROR << "matrix hasn't been encoded yet\n";
	  exit(1);
        } else { 
	  if (csx_->symmetric) {
	    if (!SetValueCsxSym<idx_t, val_t>(csx_, row_idx,
					      col_idx, value)) {
	      LOG_WARNING << "matrix entry not found\n";
	      return false;
	    }
	  } else {
	    if (!SetValueCsx<idx_t, val_t>(csx_, row_idx,
					   col_idx, value)) {
	      LOG_WARNING << "matrix entry not found\n";
	      return false;
	    }
	  }
        }

        return true;
      }

      void Save(const char *filename)
      {
        if (!csx_) {
	  LOG_ERROR << "Matrix hasn't been encoded yet\n";
	  exit(1);
        } else {
	  SaveCsx<idx_t, val_t>(csx_, filename, NULL);
        }
      }

      void Destroy()
      {
        PutSpmMt<idx_t, val_t>(csx_);
        csx_ = 0;
      }

    private:
      spm_mt_t *DoCreateCsx(internal::Sym<true>)
      {
        RtCtx& rt_context = RtCtx::GetInstance();
        SparseInternal<SparsePartitionSym<idx_t, val_t> > *spi;

        // Converting to internal representation
        timers_.StartTimer("load");
        spi = SparseInternal<SparsePartitionSym<idx_t, val_t> >::DoLoadMatrixSym
	  (*this, rt_context.GetNrThreads());
        timers_.PauseTimer("load");

        // Converting to CSX
        timers_.StartTimer("preproc");
        csx_ = BuildCsxSym<idx_t, val_t>(spi);
        timers_.PauseTimer("preproc");
        stringstream os;
        os << "\n==== GENERAL TIMING STATISTICS ====\n";
        timers_.PrintAllTimers(os);
        LOG_INFO << os.str();

        // Cleanup
        delete spi;
        return csx_;
      }

      spm_mt_t *DoCreateCsx(internal::Sym<false>)
      {
        RtCtx& rt_context = RtCtx::GetInstance();
        SparseInternal<SparsePartition<idx_t, val_t> > *spi;

        // Converting to internal representation
        timers_.StartTimer("load");
        spi = SparseInternal<SparsePartition<idx_t, val_t> >::
	  DoLoadMatrix(*this, rt_context.GetNrThreads());
        timers_.PauseTimer("load");

        // Converting to CSX
        timers_.StartTimer("preproc");
        csx_ = BuildCsx<idx_t, val_t>(spi);
        timers_.PauseTimer("preproc");
        stringstream os;
        os << "\n==== GENERAL TIMING STATISTICS ====\n";
        timers_.PrintAllTimers(os);
        LOG_VERBOSE << os.str();

        // Cleanup
        delete spi;
        return csx_;
      }
    
      void DoCreateTimers()
      {
        timers_.CreateTimer("load", "Conversion to internal representation");
        timers_.CreateTimer("preproc", "Preprocessing");
      }

      spm_mt_t *csx_;
      internal::allow_instantiation<idx_t, val_t> instance;
      timing::TimerCollection timers_;
    };

  } // end of namespace csx
} // end of namespace sparsex

#endif // SPARSEX_INTERNALS_SPARSE_MATRIX_HPP
