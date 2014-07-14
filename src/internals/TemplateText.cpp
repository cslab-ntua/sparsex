/*
 * \file TemplateText.cpp
 *
 * \brief Class for manipulating template texts
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/TemplateText.hpp>

using namespace boost;
using namespace std;

namespace sparsex {
namespace jit {

TemplateText::TemplateText(string text)
    : template_text_(text), placeholder_pattern_("\\$\\{(\\w+)\\}")
{
    ReplacePlaceholders();
}

void TemplateText::SetTemplate(string text)
{
    template_text_ = text;
    ReplacePlaceholders();
}


string TemplateText::Substitute(const std::map<string, string> &values)
{
    std::map<string, string>::const_iterator iter = values.begin();
    std::map<string, string>::const_iterator iter_end =
        values.end();

    for (; iter != iter_end; ++iter) {
        if (placeholders_.count(iter->first))
            placeholders_[iter->first] = iter->second;
        else
            LOG_ERROR << "key '" << iter->first << "' not found\n";
    }

    return DoSubstitute();
}

void TemplateText::ReplacePlaceholders()
{
    placeholders_.clear();

    // Scan template text for keys
    string::const_iterator s_start = template_text_.begin();
    string::const_iterator s_end = template_text_.end();
    match_results<string::const_iterator> match;

    while (regex_search(s_start, s_end, match, placeholder_pattern_)) {
        string key(match[1].first, match[1].second);
        placeholders_[key] = string("");
        s_start = match[0].second;
    }
}

string TemplateText::DoSubstitute()
{
    string ret;
    regex_replace(back_inserter(ret),
                  template_text_.begin(), template_text_.end(),
                  placeholder_pattern_, TextReplacer(this), match_default);
    return ret;
}

} // end of namespace jit
} // end of namespace sparsex
