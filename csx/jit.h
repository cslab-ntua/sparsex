/* -*- C++ -*-
 *
 * jit.h -- Just In Time compilation routines.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2011, Vasileios Karakasis
 * Copyright (C) 2010-2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_JIT_H__
#define CSX_JIT_H__

#include "compiler.h"
#include "CsxManager.h"
#include "encodings.h"
#include "jit.h"
#include "jit_config.h"
#include "jit_util.h"
#include "sparse_util.h"
#include "template_text.h"

#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/Module.h>
#include <llvm/LLVMContext.h>
#include <llvm/Analysis/Verifier.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/ExecutionEngine/JIT.h>
#include <llvm/LLVMContext.h>
#include <llvm/Module.h>
#include <llvm/PassManager.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Target/TargetData.h>
#include <llvm/Target/TargetOptions.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/IPO/PassManagerBuilder.h>


#include <iostream>
#include <sstream>
#include <cassert>

using namespace llvm;

namespace csx {

template<typename IndexType, typename ValueType>
class CsxManager;

/**
 *  Singleton wrapper to an LLVM execution engine.
 */
class CsxExecutionEngine {
public:
    static CsxExecutionEngine &CreateEngine()
    {
        static CsxExecutionEngine engine;
        return engine;
    }

    void AddModule(Module *mod);

    void RemoveModule(Module *mod)
    {
        llvm_engine_->removeModule(mod);
    }

    void *GetPointerToFunction(Function *fn)
    {
        return llvm_engine_->getPointerToFunction(fn);
    }

    ~CsxExecutionEngine()
    {
        if (llvm_engine_)
            delete llvm_engine_;
    };

private:
    CsxExecutionEngine() { };
    CsxExecutionEngine(const CsxExecutionEngine &); // Do not implement
    void operator=(const CsxExecutionEngine &);     // Do not implement

    //static CsxExecutionEngine *engine_;
    ExecutionEngine *llvm_engine_;
};

/// Initialize code generation. Call it once before start generating code
CsxExecutionEngine &CsxJitInit(void);

/**
 *  Code generator for a CSX matrix (one per thread)
 */
template<typename IndexType, typename ValueType>
class CsxJit {

public:
    CsxJit(CsxManager<IndexType, ValueType> *csxmg, CsxExecutionEngine *engine,
           unsigned int tid = 0, bool symmetric = false);

    ~CsxJit()
    {
        std::map<Encoding::Type, TemplateText *>::iterator it;
        for (it = mult_templates_.begin(); it != mult_templates_.end(); ++it)
            delete it->second;
        if (compiler_)
            delete compiler_;
        // if (context_)
        //     delete context_;
    }

    void GenCode(std::ostream &log);
    void *GetSpmvFn() const;

private:
    // Compile C99 source code into an LLVM module
    Module *DoCompile(const std::string &source) const;

    //
    // Obsolete -- all code optimizations are handled by ClangCompiler
    // 
    void DoOptimizeModule();

    void DoNewRowHook(std::map<std::string, std::string> &hooks,
                      std::ostream &log) const;
    void DoSpmvFnHook(std::map<std::string, std::string> &hooks,
                      std::ostream &log);

    //
    //  Code generation functions for each pattern type
    // 
    std::string DoGenDeltaCase(int delta_bits);
    std::string DoGenLinearCase(Encoding::Type type, int delta);
    std::string DoGenBlockCase(Encoding::Type type, int r, int c);
    TemplateText *GetMultTemplate(Encoding::Type type);
    
    std::string DoGenDeltaSymCase(int delta_bits);
    std::string DoGenLinearSymCase(Encoding::Type type, int delta);
    std::string DoGenBlockSymCase(Encoding::Type type, int r, int c);
    TemplateText *GetSymMultTemplate(Encoding::Type type);

    CsxManager<IndexType, ValueType> *csxmg_;
    Module *module_;
    CsxExecutionEngine *engine_;
    LLVMContext *context_;
    ClangCompiler *compiler_;
    unsigned int thread_id_;
    std::map<Encoding::Type, TemplateText *> mult_templates_;
    bool symmetric_;
    csx_t<ValueType> *csx_;
};


/*
 * Implementation of class CsxJit
 */
template<typename IndexType, typename ValueType>
CsxJit<IndexType, ValueType>::CsxJit(CsxManager<IndexType, ValueType> *csxmg,
                                     CsxExecutionEngine *engine,
                                     unsigned int tid, bool symmetric)
    : csxmg_(csxmg),
      module_(0),
      engine_(engine),
      thread_id_(tid),
      symmetric_(symmetric),
      csx_(NULL)
{
    context_ = new LLVMContext();
    compiler_ = new ClangCompiler();
    compiler_->AddHeaderSearchPath(CSX_PREFIX "/csx");
    compiler_->AddHeaderSearchPath(CSX_PREFIX "/lib/spm");
    // compiler_->SetDebugMode(true);
};

template<typename IndexType, typename ValueType>
TemplateText *CsxJit<IndexType, ValueType>::GetMultTemplate(Encoding::Type type)
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
    case Encoding::BlockRow1 ... Encoding::BlockRowMax:
        // Set the same template for all block row patterns
        mult_templates_[Encoding::BlockRow1] =
            new TemplateText(SourceFromFile(BlockRowOneTemplateSource));
        for (Encoding::Type t = Encoding::BlockRow2;
             t <= Encoding::BlockRowMax; ++t)
            mult_templates_[t] =
                new TemplateText(SourceFromFile(BlockRowTemplateSource));
        break;
    case Encoding::BlockCol1 ... Encoding::BlockColMax:
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
TemplateText *CsxJit<IndexType, ValueType>::GetSymMultTemplate(Encoding::Type type)
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
    case Encoding::BlockRow1 ... Encoding::BlockRowMax:
        // Set the same template for all block row patterns
        for (Encoding::Type t = Encoding::BlockRow1;
             t <= Encoding::BlockRowMax; ++t)
            mult_templates_[t] =
                new TemplateText(SourceFromFile(BlockRowSymTemplateSource));
        break;
    case Encoding::BlockCol1 ... Encoding::BlockColMax:
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
std::string CsxJit<IndexType, ValueType>::DoGenDeltaCase(int delta_bits)
{
    TemplateText *tmpl = GetMultTemplate(Encoding::None);
    std::map<std::string, std::string> subst_map;
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
std::string CsxJit<IndexType, ValueType>::DoGenLinearCase(Encoding::Type type,
                                                          int delta)
{
    TemplateText *tmpl = GetMultTemplate(type);
    std::map<std::string, std::string> subst_map;
    subst_map["delta"] = Stringify(delta);
    return tmpl->Substitute(subst_map);
}

template<typename IndexType, typename ValueType>
std::string CsxJit<IndexType, ValueType>::DoGenBlockCase(Encoding::Type type,
                                                         int r, int c)
{
    TemplateText *tmpl = GetMultTemplate(type);
    std::map<std::string, std::string> subst_map;
    subst_map["r"] = Stringify(r);
    subst_map["c"] = Stringify(c);
    return tmpl->Substitute(subst_map);
}

template<typename IndexType, typename ValueType>
std::string CsxJit<IndexType, ValueType>::DoGenDeltaSymCase(int delta_bits)
{
    TemplateText *tmpl = GetSymMultTemplate(Encoding::None);
    std::map<std::string, std::string> subst_map;
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
std::string CsxJit<IndexType, ValueType>::DoGenLinearSymCase(Encoding::Type type,
                                                             int delta)
{
    TemplateText *tmpl = GetSymMultTemplate(type);
    std::map<std::string, std::string> subst_map;
    subst_map["delta"] = Stringify(delta);
    return tmpl->Substitute(subst_map);
}

template<typename IndexType, typename ValueType>
std::string CsxJit<IndexType, ValueType>::DoGenBlockSymCase(Encoding::Type type,
                                                            int r, int c)
{
    TemplateText *tmpl = GetSymMultTemplate(type);
    std::map<std::string, std::string> subst_map;
    subst_map["r"] = Stringify(r);
    subst_map["c"] = Stringify(c);
    return tmpl->Substitute(subst_map);
}

template<typename IndexType, typename ValueType>
void CsxJit<IndexType, ValueType>::
DoNewRowHook(std::map<std::string, std::string> &hooks, std::ostream &log) const
{
    if (!symmetric_) {
        if (csxmg_->HasRowJmps()) {
            hooks["new_row_hook"] =
                "if (test_bit(&flags, CTL_RJMP_BIT))\n"
                "\t\t\t\ty_curr += ul_get(&ctl);\n"
                "\t\t\telse\n"
                "\t\t\t\ty_curr++;";
        } else {
            hooks["new_row_hook"] = "y_curr++;";
        }
    } else {
        if (csxmg_->HasRowJmps()) {
            hooks["new_row_hook"] =
                "if (test_bit(&flags, CTL_RJMP_BIT)) {\n"
                "\t\t\t\tint jmp = ul_get(&ctl);\n"
                "\t\t\t\tfor (i = 0; i < jmp; i++) {\n"
                "\t\t\t\t\ty[y_indx] += x[y_indx] * (*dv);\n"
                "\t\t\t\t\ty_indx++;\n"
                "\t\t\t\t\tdv++;\n"
                "\t\t\t\t}\n"
                "\t\t\t} else {\n"
                "\t\t\t\ty[y_indx] += x[y_indx] * (*dv);\n"
                "\t\t\t\ty_indx++;\n"
                "\t\t\t\tdv++;\n"
                "\t\t\t}\n";
        } else {
            hooks["new_row_hook"] =
                "y[y_indx] += x[y_indx] * (*dv);\n"
                "\t\t\ty_indx++;\n"
                "\t\t\tdv++;\n";
        }
    }
    
    if (!symmetric_) {
        if (csxmg_->HasFullColumnIndices()) {
            hooks["next_x"] = "x_curr = x + u32_get(&ctl);";
#ifdef PTR_ALIGN
            hooks["next_x"] = "align_ptr(&ctl, 4); x_curr = x + u32_get(&ctl);";
#endif
        } else {
            hooks["next_x"] = "x_curr += ul_get(&ctl);";
        }
    } else {
        if (csxmg_->HasFullColumnIndices()) {
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
DoSpmvFnHook(std::map<std::string, std::string> &hooks, std::ostream &log)
{
    std::map<long, std::string> func_entries;

    if (csxmg_) {
        typename CsxManager<IndexType, ValueType>::PatMap::iterator i_patt;
        typename CsxManager<IndexType, ValueType>::PatMap::iterator i_patt_end;
        i_patt = csxmg_->patterns.begin();
        i_patt_end = csxmg_->patterns.end();
        IndexType delta, r, c;
    
        for (; i_patt != i_patt_end; ++i_patt) {
            std::string patt_code, patt_func_entry;
            long patt_id = i_patt->second.flag;
            Encoding::Type type =
                static_cast<Encoding::Type>(i_patt->first / CSX_PID_OFFSET);
            delta = i_patt->first % CSX_PID_OFFSET;
            Encoding e(type);
            switch (e.GetType()) {
            case Encoding::None:
                assert(delta ==  8 ||
                       delta == 16 ||
                       delta == 32 ||
                       delta == 64);
                log << "type:" << e.GetFullName() << " size:" << delta
                    << " npatterns:" << i_patt->second.npatterns
                    << " nnz:" << i_patt->second.nr << std::endl;
                if (!symmetric_)
                    patt_code = DoGenDeltaCase(delta);
                else
                    patt_code = DoGenDeltaSymCase(delta);
                patt_func_entry = "delta" + Stringify(delta) + "_case";
                break;
            case Encoding::Horizontal:
                log << "type:" << e.GetFullName() << " delta:" << delta
                    << " npatterns:" << i_patt->second.npatterns
                    << " nnz:" << i_patt->second.nr << std::endl;
                if (!symmetric_)
                    patt_code = DoGenLinearCase(type, delta);
                else
                    patt_code = DoGenLinearSymCase(type, delta);
                patt_func_entry = "horiz" + Stringify(delta) + "_case";
                break;
            case Encoding::Vertical:
                log << "type:" << e.GetFullName() << " delta:" << delta
                    << " npatterns:" << i_patt->second.npatterns
                    << " nnz:" << i_patt->second.nr << std::endl;
                if (!symmetric_)
                    patt_code = DoGenLinearCase(type, delta);
                else
                    patt_code = DoGenLinearSymCase(type, delta);
                patt_func_entry = "vert" + Stringify(delta) + "_case";
                break;
            case Encoding::Diagonal:
                log << "type:" << e.GetFullName() << " delta:" << delta
                    << " npatterns:" << i_patt->second.npatterns
                    << " nnz:" << i_patt->second.nr << std::endl;
                if (!symmetric_)
                    patt_code = DoGenLinearCase(type, delta);
                else
                    patt_code = DoGenLinearSymCase(type, delta);
                patt_func_entry = "diag" + Stringify(delta) + "_case";
                break;
            case Encoding::AntiDiagonal:
                log << "type:" << e.GetFullName() << " delta:" << delta
                    << " npatterns:" << i_patt->second.npatterns
                    << " nnz:" << i_patt->second.nr << std::endl;
                if (!symmetric_)
                    patt_code = DoGenLinearCase(type, delta);
                else
                    patt_code = DoGenLinearSymCase(type, delta);
                patt_func_entry = "rdiag" + Stringify(delta) + "_case";
                break;
            case Encoding::BlockRowMin ... Encoding::BlockRowMax:
                r = e.GetBlockAlignment();
                c = delta;
                log << "type:" << e.GetFullName()
                    << " dim:" << r << "x" << c
                    << " npatterns:" << i_patt->second.npatterns
                    << " nnz:" << i_patt->second.nr << std::endl;
                if (!symmetric_)
                    patt_code = DoGenBlockCase(type, r, c);
                else
                    patt_code = DoGenBlockSymCase(type, r, c);
                patt_func_entry = "block_row_" + Stringify(r) + "x" +
                    Stringify(c) + "_case";
                break;
            case Encoding::BlockColMin ... Encoding::BlockColMax:
                r = delta;
                c = e.GetBlockAlignment();
                log << "type:" << e.GetFullName()
                    << " dim:" << r << "x" << c
                    << " npatterns:" << i_patt->second.npatterns
                    << " nnz:" << i_patt->second.nr << std::endl;
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
    } else if (csx_) {  //FIXME code duplication
#if 0
        uint8_t t;
        IndexType r, c, delta;

        for (int i = 0; csx_->id_map[i] != -1; i++) {
            t = csx_->id_map[i] / CSX_PID_OFFSET;
            delta = csx_->id_map[i] % CSX_PID_OFFSET;
            std::string patt_code, patt_func_entry;
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
            case Encoding::BlockRowMin ... Encoding::BlockRowMax:
                r = e.GetBlockAlignment();
                c = delta;
                if (!symmetric_)
                    patt_code = DoGenBlockCase(type, r, c);
                else
                    patt_code = DoGenBlockSymCase(type, r, c);
                patt_func_entry = "block_row_" + Stringify(r) + "x" +
                    Stringify(c) + "_case";
                break;
            case Encoding::BlockColMin ... Encoding::BlockColMax:
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
#endif        
	}

    // Add function table entries in the correct order
    std::map<long, std::string>::const_iterator i_fentry =
        func_entries.begin();
    std::map<long, std::string>::const_iterator fentries_end =
        func_entries.end();

    if (func_entries.size() == 1) {
        // Don't switch, just call the pattern-specific mult. routine
        if (!symmetric_) {
            hooks["body_hook"] =
                "yr += "+ Stringify(i_fentry->second) +
                "(&ctl, size, &v, &x_curr, &y_curr);";
        } else {
            hooks["body_hook"] =
                "yr += " + Stringify(i_fentry->second) +
                "(&ctl, size, &v, x, y, cur, &x_indx, &y_indx);";
        }
    } else {
        hooks["body_hook"] = "switch (patt_id) {\n";
        for (; i_fentry != fentries_end; ++i_fentry) {
            if (!symmetric_) {    
                hooks["body_hook"] +=
                    Tabify(2) + "case " + Stringify(i_fentry->first) + ":\n" +
                    Tabify(3) + "yr += " + Stringify(i_fentry->second) +
                    "(&ctl, size, &v, &x_curr, &y_curr);\n" +
                    Tabify(3) + "break;\n";
            } else {
                hooks["body_hook"] +=
                    Tabify(2) + "case " + Stringify(i_fentry->first) + ":\n" +
                    Tabify(3) + "yr += " + Stringify(i_fentry->second) +
                    "(&ctl, size, &v, x, y, cur, &x_indx, &y_indx);\n" +
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
void CsxJit<IndexType, ValueType>::GenCode(std::ostream &log)
{
    std::map<std::string, std::string> hooks;
    
    if (!symmetric_) {
        // Load the template source
        TemplateText source_tmpl(SourceFromFile(CsxTemplateSource));
        
        DoNewRowHook(hooks, log);
        DoSpmvFnHook(hooks, log);
        hooks["header_prefix"] = CSX_PREFIX;

        // Substitute and compile into an LLVM module
        compiler_->SetLogStream(&log);
        module_ = DoCompile(source_tmpl.Substitute(hooks));
        if (!module_) {
            log << "compilation failed for thread " << thread_id_ << "\n";
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
        if (!module_) {
            log << "compilation failed for thread " << thread_id_ << "\n";
            exit(1);
        }
    }
}

template<typename IndexType, typename ValueType>
Module *CsxJit<IndexType, ValueType>::DoCompile(const std::string &source) const
{
    if (compiler_->DebugMode())
        compiler_->KeepTemporaries(true);
    return compiler_->Compile(source, context_);
}

template<typename IndexType, typename ValueType>
void CsxJit<IndexType, ValueType>::DoOptimizeModule()
{
    PassManagerBuilder pm_builder;
    PassManager pm;

    pm_builder.OptLevel = 3;
    pm_builder.populateModulePassManager(pm);
    pm_builder.populateLTOPassManager(pm,
                                      true /* Internalize */,
                                      true /* RunInliner */);
    pm.run(*module_);
    return;
}

template<typename IndexType, typename ValueType>
void *CsxJit<IndexType, ValueType>::GetSpmvFn() const
{
    Function *llvm_fn;
    
    engine_->AddModule(module_);
    if (!symmetric_)
        llvm_fn = module_->getFunction("spm_csx32_double_multiply");
    else
        llvm_fn = module_->getFunction("spm_csx32_double_sym_multiply");
    assert(llvm_fn);
    return engine_->GetPointerToFunction(llvm_fn);
}

} // end csx namespace

#endif // CSX_JIT_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
