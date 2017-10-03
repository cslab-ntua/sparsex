/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Vasileios Karakasis
 * Copyright (C) 2014,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Encodings.cpp
 * \brief All about encoding types
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/Encodings.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;

namespace sparsex {
  namespace csx {

    std::map<Encoding::Type, pair<string, string> >
    Encoding::names_ = boost::assign::map_list_of
      (Encoding::None, make_pair("none", "Delta"))
      (Encoding::Horizontal, make_pair("h", "Horizontal"))
      (Encoding::Vertical, make_pair("v", "Vertical"))
      (Encoding::Diagonal, make_pair("d", "Diagonal"))
      (Encoding::AntiDiagonal, make_pair("ad", "Antidiagonal"))
      (Encoding::BlockRow1, make_pair("br1", "BlockRow1"))
      (Encoding::BlockRow2, make_pair("br2", "BlockRow2"))
      (Encoding::BlockRow3, make_pair("br3", "BlockRow3"))
      (Encoding::BlockRow4, make_pair("br4", "BlockRow4"))
      (Encoding::BlockRow5, make_pair("br5", "BlockRow5"))
      (Encoding::BlockRow6, make_pair("br6", "BlockRow6"))
      (Encoding::BlockRow7, make_pair("br7", "BlockRow7"))
      (Encoding::BlockRow8, make_pair("br8", "BlockRow8"))
      (Encoding::BlockCol1, make_pair("bc1", "BlockCol1"))
      (Encoding::BlockCol2, make_pair("bc2", "BlockCol2"))
      (Encoding::BlockCol3, make_pair("bc3", "BlockCol3"))
      (Encoding::BlockCol4, make_pair("bc4", "BlockCol4"))
      (Encoding::BlockCol5, make_pair("bc5", "BlockCol5"))
      (Encoding::BlockCol6, make_pair("bc6", "BlockCol6"))
      (Encoding::BlockCol7, make_pair("bc7", "BlockCol7"))
      (Encoding::BlockCol8, make_pair("bc8", "BlockCol8"))
      (Encoding::BlockRows, make_pair("br", "BlockRows"))
      (Encoding::BlockCols, make_pair("bc", "BlockCols"))
      (Encoding::All, make_pair("all", "all"));

    std::map<string, Encoding::Type> Encoding::InitInverseNameMap()
    {
      std::map<string, Encoding::Type> ret;

      // Index on short name
      map<Encoding::Type, pair<string, string> >::const_iterator i;
      for (i = names_.begin(); i != names_.end(); ++i) {
        Encoding::Type type = i->first;
        const string &name = (i->second).first;
        ret[name] = type;
      }

      return ret;
    }

    std::map<string, Encoding::Type> Encoding::rev_names_ =
      InitInverseNameMap();

    void Encoding::GetTypes(vector<Type> &types) const
    {
      switch (type_)
	{
	case BlockRows:
	  for (Type i = BlockRowMin; i <= BlockRowMax; ++i)
            types.push_back(i);
	  break;
	case BlockCols:
	  for (Type i = BlockColMin; i <= BlockColMax; ++i)
            types.push_back(i);
	  break;
	case All:
	  for (Type i = None; i < __EndOfTypes__; ++i)
            types.push_back(i);
	  break;
	default:
	  types.push_back(type_);
	}
    }

    void Encoding::CheckNameValidity(const string &name)
    {
      std::map<string, Type>::const_iterator i = rev_names_.find(name);
      if (i == rev_names_.end()) {
        LOG_ERROR << "invalid value \"" << name << "\" while setting property "
	  "\"spx.preproc.xform\"\n";
        exit(1);
      }
    }

    EncodingSequence::EncodingSequence(const string &seq_str)
      : explicit_(false)
    {
      regex xform_syntax("([a-z]+([0-9]*))(\\{([0-9]+(,[0-9]+)*)\\})?");
      match_results<string::const_iterator> match;
      string::const_iterator s_start = seq_str.begin();
      string::const_iterator s_end = seq_str.end();

      vector<size_t> deltas;
      while (regex_search(s_start, s_end, match, xform_syntax)) {
        string xform_str(match[1].first, match[1].second);
        // Check first if xform_str is valid
        Encoding::CheckNameValidity(xform_str);

        string deltas_str(match[4].first, match[4].second);
        s_start = match[0].second;

        // parse specific delta values
        char_separator<char> sep(",");
        tokenizer<char_separator<char> > tokens(deltas_str, sep);
        tokenizer<char_separator<char> >::iterator tok_iter;
        for (tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
	  // that's a sequence of explicit encodings
	  explicit_ = true;
	  deltas.push_back(lexical_cast<size_t>(*tok_iter));
        }

        sequence_.push_back(make_pair(Encoding(xform_str), deltas));
        deltas.clear();
      }
    }

    void EncodingSequence::Print(ostream &out) const
    {
      for (size_t i = 0; i < sequence_.size(); ++i) {
        const Encoding &enc = sequence_[i].first;
        const vector<size_t> &deltas = sequence_[i].second;
        out << enc.GetShortName();
        if (deltas.size()) {
	  out << "{";
	  for (size_t j = 0; j < deltas.size(); ++j) {
	    out << deltas[j];
	    if (j < deltas.size() - 1)
	      out << ",";
	  }

	  out << "}";
        }

        if (i < sequence_.size() - 1)
	  out << ",";
      }
    }

    bimap<PreprocMethod::Type, string> PreprocMethod::method_names_ =
      PreprocMethod::InitMethodNames();

    bimap<PreprocHeuristic::Type, string> PreprocHeuristic::heur_names_=
      PreprocHeuristic::InitHeuristicNames();

    void PreprocMethod::CheckNameValidity(const string &name)
    {
      NameMap::right_const_iterator i = method_names_.right.find(name);
      if (i == method_names_.right.end()) {
        LOG_ERROR << "invalid value \"" << name << "\" while setting property "
	  "\"spx.preproc.sampling\"\n";
        exit(1);
      }
    }

    void PreprocHeuristic::CheckNameValidity(const string &name)
    {
      NameMap::right_const_iterator i = heur_names_.right.find(name);
      if (i == heur_names_.right.end()) {
        LOG_ERROR << "invalid value \"" << name << "\" while setting property "
	  "\"spx.preproc.heuristic\"\n";
        exit(1);
      }
    }

  } // end of namespace csx
} // end of namespace sparsex
