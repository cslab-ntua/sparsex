/*
 * jit.cc -- Just In Time compilation routines.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * Copyright (C) 2011, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <iostream>
#include <sstream>
#include <cassert>

#include "llvm/Analysis/Verifier.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/Target/TargetData.h"
#include "llvm/Target/TargetOptions.h"
#include "llvm/Target/TargetSelect.h"
#include "llvm/Support/StandardPasses.h"

#include "spm.h"
#include "csx.h"
#include "drle.h"
#include "jit.h"
#include "jit_config.h"
#include "jit_util.h"

using namespace llvm;
using namespace csx;

CsxExecutionEngine &csx::CsxJitInit(void)
{
    InitializeNativeTarget();
    return CsxExecutionEngine::CreateEngine();
}

CsxJit::CsxJit(CsxManager *csxmg, CsxExecutionEngine *engine, unsigned int tid)
    : csxmg_(csxmg), module_(0), engine_(engine),
      thread_id_(tid)
{
    context_ = new LLVMContext();
    compiler_ = new ClangCompiler();
    compiler_->AddHeaderSearchPath(CSX_PREFIX "/csx");
    compiler_->AddHeaderSearchPath(CSX_PREFIX "/lib/spm");
    compiler_->AddHeaderSearchPath(CSX_PREFIX "/lib/dynarray");
};


TemplateText *CsxJit::GetMultTemplate(SpmIterOrder type)
{
    if (mult_templates_.count(type))
        goto exit;

    switch (type) {
    case NONE:
        mult_templates_[type] =
            new TemplateText(SourceFromFile(DeltaTemplateSource));
        break;
    case HORIZONTAL:
        mult_templates_[type] =
            new TemplateText(SourceFromFile(HorizTemplateSource));
        break;
    case VERTICAL:
        mult_templates_[type] =
            new TemplateText(SourceFromFile(VertTemplateSource));
        break;
    case DIAGONAL:
        mult_templates_[type] =
            new TemplateText(SourceFromFile(DiagTemplateSource));
        break;
    case REV_DIAGONAL:
        mult_templates_[type] =
            new TemplateText(SourceFromFile(RDiagTemplateSource));
        break;
    case BLOCK_R1 ... BLOCK_COL_START-1:
        // Set the same template for all block row patterns
        mult_templates_[BLOCK_R1] =
            new TemplateText(SourceFromFile(BlockRowOneTemplateSource));
        for (int t = BLOCK_R2; t < BLOCK_COL_START; ++t)
            mult_templates_[static_cast<SpmIterOrder>(t)] =
                new TemplateText(SourceFromFile(BlockRowTemplateSource));
        break;
    case BLOCK_COL_START ... BLOCK_TYPE_END-1:
        // Set the same template for all block row patterns
        mult_templates_[BLOCK_C1] =
            new TemplateText(SourceFromFile(BlockColOneTemplateSource));
        for (int t = BLOCK_C2; t < BLOCK_TYPE_END; ++t)
            mult_templates_[static_cast<SpmIterOrder>(t)] =
                new TemplateText(SourceFromFile(BlockColTemplateSource));
        break;
    default:
        assert(0 && "unknown pattern type");
    };

exit:
    return mult_templates_[type];
}

std::string CsxJit::DoGenDeltaCase(int delta_bits)
{
    TemplateText *tmpl = GetMultTemplate(NONE);
    std::map<std::string, std::string> subst_map;
    subst_map["bits"] = Stringify(delta_bits);
    int delta_bytes = delta_bits / 8;
    if (delta_bytes > 1)
        subst_map["align_ctl"] =
            "align_ptr(ctl," + Stringify(delta_bytes) + ");";

    return tmpl->Substitute(subst_map);
}

std::string CsxJit::DoGenLinearCase(SpmIterOrder type, int delta)
{
    TemplateText *tmpl = GetMultTemplate(type);
    std::map<std::string, std::string> subst_map;
    subst_map["delta"] = Stringify(delta);
    return tmpl->Substitute(subst_map);
}

std::string CsxJit::DoGenBlockCase(SpmIterOrder type, int r, int c)
{
    TemplateText *tmpl = GetMultTemplate(type);
    std::map<std::string, std::string> subst_map;
    subst_map["r"] = Stringify(r);
    subst_map["c"] = Stringify(c);
    return tmpl->Substitute(subst_map);
}

void CsxJit::DoNewRowHook(std::map<std::string, std::string> &hooks,
                          std::ostream &log) const
{
    if (csxmg_->HasRowJmps()) {
        hooks["new_row_hook"] =
            "if (test_bit(&flags, CTL_RJMP_BIT))\n"
            "\t\t\t\ty_curr += ul_get(&ctl);\n"
            "\t\t\telse\n"
            "\t\t\t\ty_curr++;";
    } else {
        hooks["new_row_hook"] = "y_curr++;";
    };

    if (csxmg_->HasFullColumnIndices()) {
        hooks["next_x"] = "x_curr = x + u32_get(&ctl);";
    } else {
        hooks["next_x"] = "x_curr += ul_get(&ctl);";
    }
}

void CsxJit::DoSpmvFnHook(std::map<std::string, std::string> &hooks,
                          std::ostream &log)
{
    CsxManager::PatMap::iterator i_patt = csxmg_->patterns.begin();
    CsxManager::PatMap::iterator i_patt_end = csxmg_->patterns.end();
    std::map<long, std::string> func_entries;

    uint64_t delta, r, c;
    for (; i_patt != i_patt_end; ++i_patt) {
        std::string patt_code, patt_func_entry;
        long patt_id = i_patt->second.flag;
        SpmIterOrder type =
            static_cast<SpmIterOrder>(i_patt->first / CSX_PID_OFFSET);
        delta = i_patt->first % CSX_PID_OFFSET;

        switch (type) {
        case NONE:
            assert(delta ==  8 ||
                   delta == 16 ||
                   delta == 32 ||
                   delta == 64);
            log << "type:DELTA size:" << delta << " npatterns:"
                << i_patt->second.npatterns << " nnz:"
                << i_patt->second.nr << std::endl;
            patt_code = DoGenDeltaCase(delta);
            patt_func_entry = "delta" + Stringify(delta) + "_case";
            break;
        case HORIZONTAL:
            log << "type:HORIZONTAL delta:" << delta
                << " npatterns:" << i_patt->second.npatterns
                << " nnz:" << i_patt->second.nr << std::endl;
            patt_code = DoGenLinearCase(type, delta);
            patt_func_entry = "horiz" + Stringify(delta) + "_case";
            break;
        case VERTICAL:
            log << "type:VERTICAL delta:" << delta
                << " npatterns:" << i_patt->second.npatterns
                << " nnz:" << i_patt->second.nr << std::endl;
            patt_code = DoGenLinearCase(type, delta);
            patt_func_entry = "vert" + Stringify(delta) + "_case";
            break;
        case DIAGONAL:
            log << "type:DIAGONAL delta:" << delta
                << " npatterns:" << i_patt->second.npatterns
                << " nnz:" << i_patt->second.nr << std::endl;
            patt_code = DoGenLinearCase(type, delta);
            patt_func_entry = "diag" + Stringify(delta) + "_case";
            break;
	case REV_DIAGONAL:
            log << "type:REV_DIAGONAL delta:" << delta
                << " npatterns:" << i_patt->second.npatterns
                << " nnz:" << i_patt->second.nr << std::endl;
            patt_code = DoGenLinearCase(type, delta);
            patt_func_entry = "rdiag" + Stringify(delta) + "_case";
            break;
        case BLOCK_R1 ... BLOCK_COL_START - 1:
            r = type - BLOCK_TYPE_START;
            c = delta;
            log << "type:" << SpmTypesNames[type]
                << " dim:" << r << "x" << c
                << " npatterns:" << i_patt->second.npatterns
                << " nnz:" << i_patt->second.nr << std::endl;
            patt_code = DoGenBlockCase(type, r, c);
            patt_func_entry = "block_row_" + Stringify(r) + "x" +
                Stringify(c) + "_case";
            break;
        case BLOCK_COL_START ... BLOCK_TYPE_END - 1:
            r = delta;
            c = type - BLOCK_COL_START;
            log << "type:" << SpmTypesNames[type] 
                << " dim:" << r << "x" << c
                << " npatterns:" << i_patt->second.npatterns
                << " nnz:" << i_patt->second.nr << std::endl;
            patt_code = DoGenBlockCase(type, r, c);
            patt_func_entry = "block_col_" + Stringify(r) + "x" +
                Stringify(c) + "_case";
            break;
        default:
            assert(0 && "unknown type");
        }

        hooks["spmv_func_definitions"] += patt_code + "\n";
        func_entries[patt_id] = patt_func_entry;
    }

    // Add function table entries in the correct order
    std::map<long, std::string>::const_iterator i_fentry =
        func_entries.begin();
    std::map<long, std::string>::const_iterator fentries_end =
        func_entries.end();

    if (func_entries.size() == 1) {
        // Don't switch, just call the pattern-specific mult. routine
        hooks["body_hook"] =
            "\t\t\tyr += " + Stringify(i_fentry->second) +
            "(&ctl, size, &v, &x_curr, &y_curr);";
    } else {
        hooks["body_hook"] = "switch (patt_id) {\n";
        for (; i_fentry != fentries_end; ++i_fentry) {
            hooks["body_hook"] +=
                "\t\tcase " + Stringify(i_fentry->first) + ":\n"
                "\t\t\tyr += " + Stringify(i_fentry->second) +
                "(&ctl, size, &v, &x_curr, &y_curr);\n"
                "\t\t\tbreak;\n";
        }
    
        // Add default statement
        hooks["body_hook"] += "\t\tdefault:\n"
            "\t\t\tfprintf(stderr, \"[BUG] unknown pattern\\n\");\n"
            "\t\t\texit(1);\n"
            "\t\t};";
    }
}

void CsxJit::GenCode(std::ostream &log)
{
    // Load the template source
    TemplateText source_tmpl(SourceFromFile(CsxTemplateSource));

    // Fill in the hooks
    std::map<std::string, std::string> hooks;
    DoNewRowHook(hooks, log);
    DoSpmvFnHook(hooks, log);

    // Substitute and compile into an LLVM module
    compiler_->SetLogStream(&log);
    //compiler_->SetDebugMode(true);
    module_ = DoCompile(source_tmpl.Substitute(hooks));
    if (!module_) {
        log << "compilation failed for thread " << thread_id_ << "\n";
        exit(1);
    }
}

Module *CsxJit::DoCompile(const std::string &source) const
{
    if (compiler_->DebugMode())
        compiler_->KeepTemporaries(true);
    return compiler_->Compile(source, context_);
}

void CsxJit::DoOptimizeModule()
{
    PassManager pm;
    pm.add(new TargetData(module_));
    createStandardModulePasses(
        &pm,
        /* -O3 */ 3,
        /* OptimizeSize */ false,
        /* UnitAtATime      */  true,
        /* UnrollLoops      */  true,
        /* SimplifyLibCalls */  true,
        /* HaveExceptions   */  false,
        /* InliningPass     */ createFunctionInliningPass());
                               
    createStandardLTOPasses(
        &pm,
        /* Internalize     */ true,
        /* RunInliner      */ true,
        /* VerifyEach      */ false);

    pm.run(*module_);
    return;
}

void *CsxJit::GetSpmvFn() const
{
    engine_->AddModule(module_);
    Function *llvm_fn = module_->getFunction("spm_csx32_double_multiply");
    assert(llvm_fn);
    return engine_->GetPointerToFunction(llvm_fn);
}

void CsxExecutionEngine::AddModule(Module *mod)
{
    if (!llvm_engine_) {
        std::string errmsg;
        llvm_engine_ = ExecutionEngine::createJIT(mod, &errmsg);
        if (!llvm_engine_) {
            std::cerr << "failed to create LLVM execution engine: "
                      << errmsg << "\n";
            exit(1);
        }
    } else {
        llvm_engine_->addModule(mod);
    }
}

#if 0
int main(int argc, char **argv)
{
    CsxManager csxmg(NULL);
    CsxJit jit(&csxmg);
    jit.GenCode(std::cout);
    return 0;
}
#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
