/*
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file TemplateText.hpp
 * \brief Class for manipulating template texts
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_TEMPLATE_TEXT_HPP
#define SPARSEX_INTERNALS_TEMPLATE_TEXT_HPP

#include <sparsex/internals/logger/Logger.hpp>
#include <map>
#include <string>
#include <boost/regex.hpp>

namespace sparsex {
  namespace jit {

    class TemplateText
    {
    public:
      TemplateText(string text);

      const string &GetTemplate() const
      {
        return const_cast<string&>(template_text_);
      }

      void SetTemplate(string text);
      string Substitute(const std::map<string, string> &values);

    private:
      void ReplacePlaceholders();
      string DoSubstitute();
      string template_text_;
      std::map<string, string> placeholders_;
      boost::regex placeholder_pattern_;

      class TextReplacer {
      public:
        TextReplacer(TemplateText *outer)
	  : outer_(outer) {}

        string operator()(const boost::match_results<string::iterator> &match)
        {
	  string key(match[1].first, match[1].second);
	  return outer_->placeholders_[key];
        }

      private:
        TemplateText *outer_;
      };
    };

  } // end of namespace jit
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_TEMPLATE_TEXT_HPP
