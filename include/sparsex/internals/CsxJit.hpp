/*
 * Copyright (C) 2009-2017, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2014, Vasileios Karakasis
 * Copyright (C) 2010-2011, Theodoros Gkountouvas
 * Copyright (C) 2014-2017, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CsxJit.hpp
 * \brief Just In Time compilation of CSX SpMV routines
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2017
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_JIT_HPP
#define SPARSEX_INTERNALS_CSX_JIT_HPP

#include <sparsex/internals/CodeExecutor.hpp>
#include <sparsex/internals/CsxManager.hpp>
#include <sparsex/internals/Encodings.hpp>
#include <sparsex/internals/Element.hpp>
#include <sparsex/internals/JitCompiler.hpp>
#include <sparsex/internals/JitConfig.hpp>
#include <sparsex/internals/JitUtil.hpp>
#include <sparsex/internals/SpmvMethod.hpp>
#include <sparsex/internals/TemplateText.hpp>

#include <llvm/Support/raw_os_ostream.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/IR/Verifier.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <cassert>

using namespace std;
using namespace sparsex::csx;

namespace sparsex {

  // Forward declaration
  namespace csx {
    template<typename IndexType, typename ValueType>
    class CsxManager;
  }

  namespace jit {

    /**
     *  Code generator for a CSX matrix (one per thread)
     */
    template<typename IndexType, typename ValueType>
    class CsxJit {

    public:
      CsxJit(CsxManager<IndexType, ValueType> *mng,
	     unsigned int tid = 0, bool symmetric = false);

      CsxJit(CsxMatrix<IndexType, ValueType> *csx,
	     unsigned int tid, bool symmetric, 
	     bool full_column_indices, bool row_jumps);

      ~CsxJit()
      {
        std::map<Encoding::Type, TemplateText *>::iterator it;
        for (it = mult_templates_.begin(); it != mult_templates_.end(); ++it)
	  delete it->second;
      }

      void GenCode(ostream &log);
      spmv_fn_t GetSpmvFn();

    private:
      // Compile C99 source code into an LLVM module
      unique_ptr<llvm::Module> DoCompile(const string &source) const;

      void DoNewRowHook(std::map<string, string> &hooks, ostream &log) const;
      void DoSpmvFnHook(std::map<string, string> &hooks, ostream &log);

      //
      //  Code generation functions for each pattern type
      // 
      string DoGenDeltaCase(int delta_bits);
      string DoGenLinearCase(Encoding::Type type, int delta);
      string DoGenBlockCase(Encoding::Type type, int r, int c);
      TemplateText *GetMultTemplate(Encoding::Type type);
    
      string DoGenDeltaSymCase(int delta_bits);
      string DoGenLinearSymCase(Encoding::Type type, int delta);
      string DoGenBlockSymCase(Encoding::Type type, int r, int c);
      TemplateText *GetSymMultTemplate(Encoding::Type type);

      CsxManager<IndexType, ValueType> *mng_;
      unique_ptr<ClangCompiler> compiler_;
      unique_ptr<llvm::Module> module_;
      unsigned int thread_id_;
      std::map<Encoding::Type, TemplateText *> mult_templates_;
      bool symmetric_, full_column_indices_, row_jumps_;
      CsxMatrix<IndexType, ValueType> *csx_;
    };

    /*
     * Implementation of class CsxJit
     */
    template<typename IndexType, typename ValueType>
    CsxJit<IndexType, ValueType>::
    CsxJit(CsxManager<IndexType, ValueType> *mng, unsigned int tid,
	   bool symmetric)
      : mng_(mng),
	module_(),
	thread_id_(tid),
	symmetric_(symmetric),
	full_column_indices_(mng->HasFullColumnIndices()),
	row_jumps_(mng->HasRowJmps()),
	csx_(0)
    {
      compiler_.reset(new ClangCompiler());
      // FIXME
      // compiler_->AddIncludeSearchPath(SPX_JIT_INCLUDE,
      //                                 ClangCompiler::IncludePathUser);
      // compiler_->SetDebugMode(true);
    }

    template<typename IndexType, typename ValueType>
    CsxJit<IndexType, ValueType>::
    CsxJit(CsxMatrix<IndexType, ValueType> *csx, unsigned int tid,
	   bool symmetric, bool full_column_indices, bool row_jumps)
      : mng_(0),
	module_(),
	thread_id_(tid),
	symmetric_(symmetric),
	full_column_indices_(full_column_indices),
	row_jumps_(row_jumps),
	csx_(csx)
    {
      compiler_.reset(new ClangCompiler());
      // compiler_->AddIncludeSearchPath(SPX_JIT_INCLUDE,
      //                                 ClangCompiler::IncludePathUser);
      // compiler_->SetDebugMode(true);
    }

    template<typename IndexType, typename ValueType>
    TemplateText *CsxJit<IndexType, ValueType>::
    GetMultTemplate(Encoding::Type type)
    {
      if (mult_templates_.count(type))
        goto exit;

      switch (type) {
      case Encoding::None:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(DeltaTemplateSource));
        break;
      case Encoding::Horizontal:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(HorizTemplateSource));
        break;
      case Encoding::Vertical:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(VertTemplateSource));
        break;
      case Encoding::Diagonal:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(DiagTemplateSource));
        break;
      case Encoding::AntiDiagonal:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(RDiagTemplateSource));
        break;
      case Encoding::BlockRow1:
      case Encoding::BlockRow2:
      case Encoding::BlockRow3:
      case Encoding::BlockRow4:
      case Encoding::BlockRow5:
      case Encoding::BlockRow6:
      case Encoding::BlockRow7:
      case Encoding::BlockRow8:
        // Set the same template for all block row patterns
        mult_templates_[Encoding::BlockRow1] =
	  new TemplateText(SourceFromFile(BlockRowOneTemplateSource));
        for (Encoding::Type t = Encoding::BlockRow2;
             t <= Encoding::BlockRowMax; ++t)
	  mult_templates_[t] =
	    new TemplateText(SourceFromFile(BlockRowTemplateSource));
        break;
      case Encoding::BlockCol1:
      case Encoding::BlockCol2:
      case Encoding::BlockCol3:
      case Encoding::BlockCol4:
      case Encoding::BlockCol5:
      case Encoding::BlockCol6:
      case Encoding::BlockCol7:
      case Encoding::BlockCol8:
        // Set the same template for all block row patterns
        mult_templates_[Encoding::BlockCol1] =
	  new TemplateText(SourceFromFile(BlockColOneTemplateSource));
        for (Encoding::Type t = Encoding::BlockCol2;
             t <= Encoding::BlockColMax; ++t)
	  mult_templates_[t] =
	    new TemplateText(SourceFromFile(BlockColTemplateSource));
        break;
      default:
        assert(0 && "unknown pattern type");
      };

    exit:
      return mult_templates_[type];
    }

    template<typename IndexType, typename ValueType>
    TemplateText *CsxJit<IndexType, ValueType>::
    GetSymMultTemplate(Encoding::Type type)
    {
      if (mult_templates_.count(type))
        goto exit;

      switch (type) {
      case Encoding::None:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(DeltaSymTemplateSource));
        break;
      case Encoding::Horizontal:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(HorizSymTemplateSource));
        break;
      case Encoding::Vertical:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(VertSymTemplateSource));
        break;
      case Encoding::Diagonal:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(DiagSymTemplateSource));
        break;
      case Encoding::AntiDiagonal:
        mult_templates_[type] =
	  new TemplateText(SourceFromFile(RDiagSymTemplateSource));
        break;
      case Encoding::BlockRow1:
      case Encoding::BlockRow2:
      case Encoding::BlockRow3:
      case Encoding::BlockRow4:
      case Encoding::BlockRow5:
      case Encoding::BlockRow6:
      case Encoding::BlockRow7:
      case Encoding::BlockRow8:
        // Set the same template for all block row patterns
        for (Encoding::Type t = Encoding::BlockRow1;
             t <= Encoding::BlockRowMax; ++t)
	  mult_templates_[t] =
	    new TemplateText(SourceFromFile(BlockRowSymTemplateSource));
        break;
      case Encoding::BlockCol1:
      case Encoding::BlockCol2:
      case Encoding::BlockCol3:
      case Encoding::BlockCol4:
      case Encoding::BlockCol5:
      case Encoding::BlockCol6:
      case Encoding::BlockCol7:
      case Encoding::BlockCol8:
        // Set the same template for all block row patterns
        for (Encoding::Type t = Encoding::BlockCol1;
             t <= Encoding::BlockColMax; ++t)
	  mult_templates_[t] =
	    new TemplateText(SourceFromFile(BlockColSymTemplateSource));
        break;
      default:
        assert(0 && "unknown pattern type");
      };

    exit:
      return mult_templates_[type];
    }

    template<typename IndexType, typename ValueType>
    string CsxJit<IndexType, ValueType>::
    DoGenDeltaCase(int delta_bits)
    {
      TemplateText *tmpl = GetMultTemplate(Encoding::None);
      std::map<string, string> subst_map;
      subst_map["bits"] = Stringify(delta_bits);
#ifdef PTR_ALIGN
      int delta_bytes = delta_bits / 8;
      if (delta_bytes > 1)
        subst_map["align_ctl"] = 
	  "align_ptr(ctl," + Stringify(delta_bytes) + ");";
      else
        subst_map["align_ctl"] = "";
#endif
      return tmpl->Substitute(subst_map);
    }

    template<typename IndexType, typename ValueType>
    string CsxJit<IndexType, ValueType>::
    DoGenLinearCase(Encoding::Type type, int delta)
    {
      TemplateText *tmpl = GetMultTemplate(type);
      std::map<string, string> subst_map;
      subst_map["delta"] = Stringify(delta);
      return tmpl->Substitute(subst_map);
    }

    template<typename IndexType, typename ValueType>
    string CsxJit<IndexType, ValueType>::
    DoGenBlockCase(Encoding::Type type,	int r, int c)
    {
      TemplateText *tmpl = GetMultTemplate(type);
      std::map<string, string> subst_map;
      subst_map["r"] = Stringify(r);
      subst_map["c"] = Stringify(c);
      return tmpl->Substitute(subst_map);
    }

    template<typename IndexType, typename ValueType>
    string CsxJit<IndexType, ValueType>::
    DoGenDeltaSymCase(int delta_bits)
    {
      TemplateText *tmpl = GetSymMultTemplate(Encoding::None);
      std::map<string, string> subst_map;
      subst_map["bits"] = Stringify(delta_bits);
#ifdef PTR_ALIGN
      int delta_bytes = delta_bits / 8;
      if (delta_bytes > 1)
        subst_map["align_ctl"] =
	  "align_ptr(ctl," + Stringify(delta_bytes) + ");";
      else
        subst_map["align_ctl"] = "";
#endif
      return tmpl->Substitute(subst_map);
    }

    template<typename IndexType, typename ValueType>
    string CsxJit<IndexType, ValueType>::
    DoGenLinearSymCase(Encoding::Type type, int delta)
    {
      TemplateText *tmpl = GetSymMultTemplate(type);
      std::map<string, string> subst_map;
      subst_map["delta"] = Stringify(delta);
      return tmpl->Substitute(subst_map);
    }

    template<typename IndexType, typename ValueType>
    string CsxJit<IndexType, ValueType>::
    DoGenBlockSymCase(Encoding::Type type, int r, int c)
    {
      TemplateText *tmpl = GetSymMultTemplate(type);
      std::map<string, string> subst_map;
      subst_map["r"] = Stringify(r);
      subst_map["c"] = Stringify(c);
      return tmpl->Substitute(subst_map);
    }

    template<typename IndexType, typename ValueType>
    void CsxJit<IndexType, ValueType>::
    DoNewRowHook(std::map<string, string> &hooks, ostream &log) const
    {
      if (!symmetric_) {
        if (row_jumps_) {
	  hooks["new_row_hook"] =
	    "if (test_bit(&flags, CTL_RJMP_BIT))\n"
	    "\t\t\t\ty_curr += ul_get(&ctl);\n"
	    "\t\t\telse\n"
	    "\t\t\t\ty_curr++;";
        } else {
	  hooks["new_row_hook"] = "y_curr++;";
        }
      } else {
        if (row_jumps_) {
	  hooks["new_row_hook"] =
	    "if (test_bit(&flags, CTL_RJMP_BIT)) {\n"
	    "\t\t\t\tint jmp = ul_get(&ctl);\n"
	    "\t\t\t\tfor (i = 0; i < jmp; i++) {\n"
	    "\t\t\t\t\ty[y_indx] += x[y_indx] * (*dv) * scale_f;\n"
	    "\t\t\t\t\ty_indx++;\n"
	    "\t\t\t\t\tdv++;\n"
	    "\t\t\t\t}\n"
	    "\t\t\t} else {\n"
	    "\t\t\t\ty[y_indx] += x[y_indx] * (*dv) * scale_f;\n"
	    "\t\t\t\ty_indx++;\n"
	    "\t\t\t\tdv++;\n"
	    "\t\t\t}\n";
        } else {
	  hooks["new_row_hook"] =
	    "y[y_indx] += x[y_indx] * (*dv) * scale_f;\n"
	    "\t\t\ty_indx++;\n"
	    "\t\t\tdv++;\n";
        }
      }
    
      if (!symmetric_) {
        if (full_column_indices_) {
	  hooks["next_x"] = "x_curr = x + u32_get(&ctl);";
#ifdef PTR_ALIGN
	  hooks["next_x"] = "align_ptr(&ctl, 4); x_curr = x + u32_get(&ctl);";
#endif
        } else {
	  hooks["next_x"] = "x_curr += ul_get(&ctl);";
        }
      } else {
        if (full_column_indices_) {
	  hooks["next_x"] = "x_indx = u32_get(&ctl);";
#if PTR_ALIGN
	  hooks["next_x"] = "align_ptr(&ctl, 4); x_curr = x + u32_get(&ctl);";
#endif
        } else {
	  hooks["next_x"] = "x_indx += ul_get(&ctl);";
        }
      }
    }

    template<typename IndexType, typename ValueType>
    void CsxJit<IndexType, ValueType>::
    DoSpmvFnHook(std::map<string, string> &hooks, ostream &log)
    {
      std::map<long, string> func_entries;

      if (mng_) {
        typename CsxManager<IndexType, ValueType>::PatMap::iterator i_patt;
        typename CsxManager<IndexType, ValueType>::PatMap::iterator i_patt_end;
        i_patt = mng_->patterns.begin();
        i_patt_end = mng_->patterns.end();
        IndexType delta, r, c;
    
        for (; i_patt != i_patt_end; ++i_patt) {
	  string patt_code, patt_func_entry;
	  long patt_id = i_patt->second.flag;
	  Encoding::Type type =
	    static_cast<Encoding::Type>(i_patt->first / PatternIdOffset);
	  delta = i_patt->first % PatternIdOffset;
	  Encoding e(type);
	  switch (e.GetType()) {
	  case Encoding::None:
	    assert(delta ==  8 ||
		   delta == 16 ||
		   delta == 32 ||
		   delta == 64);
	    log << "type:" << e.GetFullName() << " size:" << delta
		<< " npatterns:" << i_patt->second.npatterns
		<< " nnz:" << i_patt->second.nr << endl;
	    if (!symmetric_)
	      patt_code = DoGenDeltaCase(delta);
	    else
	      patt_code = DoGenDeltaSymCase(delta);
	    patt_func_entry = "delta" + Stringify(delta) + "_case";
	    break;
	  case Encoding::Horizontal:
	    log << "type:" << e.GetFullName() << " delta:" << delta
		<< " npatterns:" << i_patt->second.npatterns
		<< " nnz:" << i_patt->second.nr << endl;
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "horiz" + Stringify(delta) + "_case";
	    break;
	  case Encoding::Vertical:
	    log << "type:" << e.GetFullName() << " delta:" << delta
		<< " npatterns:" << i_patt->second.npatterns
		<< " nnz:" << i_patt->second.nr << endl;
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "vert" + Stringify(delta) + "_case";
	    break;
	  case Encoding::Diagonal:
	    log << "type:" << e.GetFullName() << " delta:" << delta
		<< " npatterns:" << i_patt->second.npatterns
		<< " nnz:" << i_patt->second.nr << endl;
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "diag" + Stringify(delta) + "_case";
	    break;
	  case Encoding::AntiDiagonal:
	    log << "type:" << e.GetFullName() << " delta:" << delta
		<< " npatterns:" << i_patt->second.npatterns
		<< " nnz:" << i_patt->second.nr << endl;
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "rdiag" + Stringify(delta) + "_case";
	    break;
	  case Encoding::BlockRow1:
	  case Encoding::BlockRow2:
	  case Encoding::BlockRow3:
	  case Encoding::BlockRow4:
	  case Encoding::BlockRow5:
	  case Encoding::BlockRow6:
	  case Encoding::BlockRow7:
	  case Encoding::BlockRow8:
	    r = e.GetBlockAlignment();
	    c = delta;
	    log << "type:" << e.GetFullName()
		<< " dim:" << r << "x" << c
		<< " npatterns:" << i_patt->second.npatterns
		<< " nnz:" << i_patt->second.nr << endl;
	    if (!symmetric_)
	      patt_code = DoGenBlockCase(type, r, c);
	    else
	      patt_code = DoGenBlockSymCase(type, r, c);
	    patt_func_entry = "block_row_" + Stringify(r) + "x" +
	      Stringify(c) + "_case";
	    break;
	  case Encoding::BlockCol1:
	  case Encoding::BlockCol2:
	  case Encoding::BlockCol3:
	  case Encoding::BlockCol4:
	  case Encoding::BlockCol5:
	  case Encoding::BlockCol6:
	  case Encoding::BlockCol7:
	  case Encoding::BlockCol8:
	    r = delta;
	    c = e.GetBlockAlignment();
	    log << "type:" << e.GetFullName()
		<< " dim:" << r << "x" << c
		<< " npatterns:" << i_patt->second.npatterns
		<< " nnz:" << i_patt->second.nr << endl;
	    if (!symmetric_)
	      patt_code = DoGenBlockCase(type, r, c);
	    else
	      patt_code = DoGenBlockSymCase(type, r, c);
	    patt_func_entry = "block_col_" + Stringify(r) + "x" +
	      Stringify(c) + "_case";
	    break;
	  default:
	    assert(0 && "unknown type");
	  }

	  hooks["spmv_func_definitions"] += patt_code + "\n";
	  func_entries[patt_id] = patt_func_entry;
        }
      } else if (csx_) {
        IndexType r, c, delta;

        for (int i = 0; csx_->id_map[i] != -1; i++) {
	  Encoding::Type type =
	    static_cast<Encoding::Type>(csx_->id_map[i] / PatternIdOffset);
	  delta = csx_->id_map[i] % PatternIdOffset;
	  string patt_code, patt_func_entry;
	  Encoding e(type);
	  switch (e.GetType()) {
	  case Encoding::None:
	    assert(delta ==  8 ||
		   delta == 16 ||
		   delta == 32 ||
		   delta == 64);
	    if (!symmetric_)
	      patt_code = DoGenDeltaCase(delta);
	    else
	      patt_code = DoGenDeltaSymCase(delta);
	    patt_func_entry = "delta" + Stringify(delta) + "_case";
	    break;
	  case Encoding::Horizontal:
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "horiz" + Stringify(delta) + "_case";
	    break;
	  case Encoding::Vertical:
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "vert" + Stringify(delta) + "_case";
	    break;
	  case Encoding::Diagonal:
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "diag" + Stringify(delta) + "_case";
	    break;
	  case Encoding::AntiDiagonal:
	    if (!symmetric_)
	      patt_code = DoGenLinearCase(type, delta);
	    else
	      patt_code = DoGenLinearSymCase(type, delta);
	    patt_func_entry = "rdiag" + Stringify(delta) + "_case";
	    break;
	  case Encoding::BlockRow1:
	  case Encoding::BlockRow2:
	  case Encoding::BlockRow3:
	  case Encoding::BlockRow4:
	  case Encoding::BlockRow5:
	  case Encoding::BlockRow6:
	  case Encoding::BlockRow7:
	  case Encoding::BlockRow8:
	    r = e.GetBlockAlignment();
	    c = delta;
	    if (!symmetric_)
	      patt_code = DoGenBlockCase(type, r, c);
	    else
	      patt_code = DoGenBlockSymCase(type, r, c);
	    patt_func_entry = "block_row_" + Stringify(r) + "x" +
	      Stringify(c) + "_case";
	    break;
	  case Encoding::BlockCol1:
	  case Encoding::BlockCol2:
	  case Encoding::BlockCol3:
	  case Encoding::BlockCol4:
	  case Encoding::BlockCol5:
	  case Encoding::BlockCol6:
	  case Encoding::BlockCol7:
	  case Encoding::BlockCol8:
	    r = delta;
	    c = e.GetBlockAlignment();
	    if (!symmetric_)
	      patt_code = DoGenBlockCase(type, r, c);
	    else
	      patt_code = DoGenBlockSymCase(type, r, c);
	    patt_func_entry = "block_col_" + Stringify(r) + "x" +
	      Stringify(c) + "_case";
	    break;
	  default:
	    assert(0 && "unknown type");
	  }

	  hooks["spmv_func_definitions"] += patt_code + "\n";
	  func_entries[i] = patt_func_entry;
        }
      }           

      // Add function table entries in the correct order
      std::map<long, string>::const_iterator i_fentry = func_entries.begin();
      std::map<long, string>::const_iterator fentries_end = func_entries.end();

      if (func_entries.size() == 1) {
        // Don't switch, just call the pattern-specific mult. routine
        if (!symmetric_) {
	  hooks["body_hook"] =
	    "yr += "+ Stringify(i_fentry->second) +
	    "(&ctl, size, &v, &x_curr, &y_curr, scale_f);";
        } else {
	  hooks["body_hook"] =
	    "yr += " + Stringify(i_fentry->second) +
	    "(&ctl, size, &v, x, y, cur, &x_indx, &y_indx, scale_f);";
        }
      } else {
        hooks["body_hook"] = "switch (patt_id) {\n";
        for (; i_fentry != fentries_end; ++i_fentry) {
	  if (!symmetric_) {    
	    hooks["body_hook"] +=
	      Tabify(2) + "case " + Stringify(i_fentry->first) + ":\n" +
	      Tabify(3) + "yr += " + Stringify(i_fentry->second) +
	      "(&ctl, size, &v, &x_curr, &y_curr, scale_f);\n" +
	      Tabify(3) + "break;\n";
	  } else {
	    hooks["body_hook"] +=
	      Tabify(2) + "case " + Stringify(i_fentry->first) + ":\n" +
	      Tabify(3) + "yr += " + Stringify(i_fentry->second) +
	      "(&ctl, size, &v, x, y, cur, &x_indx, &y_indx, scale_f);\n" +
	      Tabify(3) + "break;\n";
	  }
        }
        
        // Add default statement
        hooks["body_hook"] +=
	  Tabify(2) + "default:\n" +
	  Tabify(3) + "fprintf(stderr, \"[BUG] unknown pattern\\n\");\n" +
	  Tabify(3) + "exit(1);\n" +
	  Tabify(2) + "};";
      }
    }

    template<typename IndexType, typename ValueType>
    void CsxJit<IndexType, ValueType>::
    GenCode(ostream &log)
    {
      std::map<string, string> hooks;
    
      if (!symmetric_) {
        // Load the template source
        TemplateText source_tmpl(SourceFromFile(CsxTemplateSource));
        
        DoNewRowHook(hooks, log);
        DoSpmvFnHook(hooks, log);

        // Substitute and compile into an LLVM module
        compiler_->SetLogStream(&log);
        module_ = DoCompile(source_tmpl.Substitute(hooks));
	if (module_.get() == nullptr) {
	  log << "[ERROR]: dynamic compilation failed for thread "
	      << thread_id_ << "\n";
	  exit(1);
	}
      } else {
        // Load the template source
        TemplateText source_tmpl(SourceFromFile(CsxSymTemplateSource));    
        DoNewRowHook(hooks, log);
        DoSpmvFnHook(hooks, log);

        // Substitute and compile into an LLVM module
        compiler_->SetLogStream(&log);
        module_ = DoCompile(source_tmpl.Substitute(hooks));
	if (module_.get() == nullptr) {
	  log << "[ERROR]: dynamic compilation failed for thread "
	      << thread_id_ << "\n";
	  exit(1);
	}
      }

      // FIXME
      // // Check whether the compiled function is in the current module and return
      // // its mangled name
      // if (!GetMangledName(info_.mangled_name, module_.get(), info_.func_name)) {
      //   cerr << "[ERROR]: unable to find mangled name of function '"
      // 	   << info_.func_name
      // 	   << "' in the LLVM module"
      // 	   << endl;
      //   exit(1);
      // }

      // // Validate the generated code, checking for consistency
      // llvm::raw_os_ostream os(cerr);
      // if (llvm::verifyFunction(*module_->getFunction(info_.mangled_name), &os)) {
      //   cerr << "[ERROR]: failed to verify function '"
      // 	   << info_.func_name
      // 	   << "' in the LLVM module"
      // 	   << endl;
      //   exit(1);
      // }
    }

    template<typename IndexType, typename ValueType>
    unique_ptr<llvm::Module> CsxJit<IndexType, ValueType>::
    DoCompile(const string &source) const
    {
      compiler_->SetLogStream(&cerr);
      if (compiler_->DebugMode())
	compiler_->KeepTemporaries(true);
      return compiler_->Compile(source);
    }
  
    template<typename IndexType, typename ValueType>
    spmv_fn_t CsxJit<IndexType, ValueType>::
    GetSpmvFn()
    { 
      CodeExecutor &executor = CodeExecutor::GetInstance();
      executor.AddModule(move(module_));
      string llvm_fn;
      if (!symmetric_)
	llvm_fn = "spm_csx_multiply";
      else
	llvm_fn = "spm_csx_sym_multiply";
      spmv_fn_t fn_ptr =
	reinterpret_cast<spmv_fn_t>(executor.GetFnAddress(llvm_fn));
      if (!fn_ptr) {
	cerr << "[ERROR]: unable to find function '"
	     << llvm_fn
	     << "' in execution engine"
	     << endl;
	exit(1);
      }

      return fn_ptr;
    }

  } // end of namespace jit
} // end of namespace sparsex

#endif // SPARSEX_INTERNALS_CSX_JIT_HPP
