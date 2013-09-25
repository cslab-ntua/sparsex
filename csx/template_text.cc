/*
 * template_text.cc -- Class for manipulating template texts.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "template_text.h"
#include <boost/regex.hpp>

using namespace boost;

TemplateText::TemplateText(std::string text)
    : template_text_(text), placeholder_pattern_("\\$\\{(\\w+)\\}")
{
    ReplacePlaceholders();
}

void TemplateText::SetTemplate(std::string text)
{
    template_text_ = text;
    ReplacePlaceholders();
}


std::string TemplateText::Substitute(const std::map<std::string, std::string>
                                     &values)
{
    std::map<std::string, std::string>::const_iterator iter = values.begin();
    std::map<std::string, std::string>::const_iterator iter_end =
        values.end();

    for (; iter != iter_end; ++iter) {
        if (placeholders_.count(iter->first))
            placeholders_[iter->first] = iter->second;
        else
            std::cerr << "key '" << iter->first << "' not found\n";
    }

    return DoSubstitute();
}

void TemplateText::ReplacePlaceholders()
{
    placeholders_.clear();

    // Scan template text for keys
    std::string::const_iterator s_start = template_text_.begin();
    std::string::const_iterator s_end = template_text_.end();
    match_results<std::string::const_iterator> match;

    while (regex_search(s_start, s_end, match, placeholder_pattern_)) {
        std::string key(match[1].first, match[1].second);
        placeholders_[key] = std::string("");
        s_start = match[0].second;
    }
}

std::string TemplateText::DoSubstitute()
{
    std::string ret;
    regex_replace(back_inserter(ret),
                  template_text_.begin(), template_text_.end(),
                  placeholder_pattern_, TextReplacer(this), match_default);
    return ret;
}
