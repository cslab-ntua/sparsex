/*
 * Copyright (C) 2012-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2014, Athena Elafrou
 * Copyright (C) 2013,      Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Runtime.cpp
 * \brief Front-end utilities for runtime configuration
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Runtime.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/assign/list_of.hpp>

using namespace std;
using namespace boost;

namespace sparsex {
  namespace runtime {

    /**
     *  Load default CSX runtime properties
     */
    RtConfig::PropertyMap RtConfig::DefaultProperties()
    {
      return boost::assign::map_list_of
        (RtNrThreads, "1")
        (RtCpuAffinity, "0")
#if SPX_USE_NUMA
        (PreprocHeuristic, "cost")
#else
        (PreprocHeuristic, "ratio")
#endif
        (PreprocXform, "all")
        (PreprocMethod, "portion")
        (PreprocNrSamples, "48")
        (PreprocSamplingPortion, "0.01")
        (PreprocWindowSize, "0")
        (MatrixSymmetric, "false")
        (MatrixSplitBlocks, "true")
        (MatrixOneDimBlocks, "false")
#if SPX_USE_NUMA
        (MatrixFullColind, "true")
#else
        (MatrixFullColind, "false")
#endif            
        (MatrixMinUnitSize, "4")
        (MatrixMaxUnitSize, "255")
        (MatrixMinCoverage, "0.1");
    }

    void RtConfig::InitPropertyMnemonics()
    {
      mnemonic_map_.insert(mnemonic_pair(RtNrThreads, 
					 "spx.rt.nr_threads"));
      mnemonic_map_.insert(mnemonic_pair(RtCpuAffinity, 
					 "spx.rt.cpu_affinity"));
      mnemonic_map_.insert(mnemonic_pair(PreprocHeuristic, 
					 "spx.preproc.heuristic"));
      mnemonic_map_.insert(mnemonic_pair(PreprocXform, 
					 "spx.preproc.xform"));
      mnemonic_map_.insert(mnemonic_pair(PreprocMethod, 
					 "spx.preproc.sampling"));
      mnemonic_map_.insert(mnemonic_pair(PreprocNrSamples,
					 "spx.preproc.sampling.nr_samples"));
      mnemonic_map_.insert(mnemonic_pair(PreprocSamplingPortion,
					 "spx.preproc.sampling.portion"));
      mnemonic_map_.insert(mnemonic_pair(PreprocWindowSize,
					 "spx.preproc.sampling.window_size"));
      mnemonic_map_.insert(mnemonic_pair(MatrixSymmetric,
					 "spx.matrix.symmetric"));
      mnemonic_map_.insert(mnemonic_pair(MatrixSplitBlocks, 
					 "spx.matrix.split_blocks"));
      mnemonic_map_.insert(mnemonic_pair(MatrixFullColind,
					 "spx.matrix.full_colind"));
      mnemonic_map_.insert(mnemonic_pair(MatrixMinUnitSize,
					 "spx.matrix.min_unit_size"));
      mnemonic_map_.insert(mnemonic_pair(MatrixMaxUnitSize,
					 "spx.matrix.max_unit_size"));
      mnemonic_map_.insert(mnemonic_pair(MatrixMinCoverage,
					 "spx.matrix.min_coverage"));
    }

    RtConfig &RtConfig::LoadFromEnv()
    {
      const char *mat_str = getenv("SYMMETRIC");
      if (mat_str) {
        SetProperty(RtConfig::MatrixSymmetric,
                    string(mat_str));
      }

      const char *mt_conf_str = getenv("CPU_AFFINITY");
      if (mt_conf_str) {
        SetProperty(RtConfig::RtCpuAffinity, string(mt_conf_str));
      }

      const char *thr_str = getenv("NUM_THREADS");
      if (thr_str) {
        SetProperty(RtConfig::RtNrThreads, string(thr_str));
      }

      const char *xform_conf_str = getenv("XFORM_CONF");
      if (xform_conf_str) {
        SetProperty(RtConfig::PreprocXform, string(xform_conf_str));
      }

      const char *wsize_str = getenv("WINDOW_SIZE");
      if (wsize_str) {
        // Automatically enable sampling
        SetProperty(RtConfig::PreprocMethod, "window");
        SetProperty(RtConfig::PreprocWindowSize, string(wsize_str));
      }

      const char *samples_str = getenv("SAMPLES");
      if (samples_str) {
        // Automatically enable sampling
        SetProperty(RtConfig::PreprocNrSamples,
                    string(samples_str));
      }

      const char *sampling_portion_str = getenv("SAMPLING_PORTION");
      if (sampling_portion_str) {
        // Automatically enable sampling
        SetProperty(RtConfig::PreprocMethod, "portion");
        SetProperty(RtConfig::PreprocSamplingPortion,
                    string(sampling_portion_str));
      }

      const char *sampling_str = getenv("SAMPLING");
      if (sampling_str) {
        SetProperty(RtConfig::PreprocMethod, 
                    string(sampling_str));
      }

      return *this;
    }

    vector<size_t> &ParseCpuAffinity(string str, vector<size_t> &affinity)
    {
      vector<string> mt_split;
      boost::split(mt_split, str, boost::algorithm::is_any_of(","),
		   boost::algorithm::token_compress_on);

      // Clear in case it has been already set.
      affinity.clear();
      for (size_t i = 0; i < mt_split.size(); ++i) {
        try {
	  affinity.push_back(boost::lexical_cast<size_t, string>
			     (mt_split[i]));
        } catch (const boost::bad_lexical_cast &e) {
	  LOG_ERROR << "invalid value \"" << mt_split[i] 
		    << "\" while setting property \"spx.rt.cpu_affinity\"\n";
	  exit(1);
        }
      }
    
      return affinity;
    }

    void RtCtx::CheckParams(const RtConfig &conf)
    {
      vector<size_t> affinity;
      affinity = ParseCpuAffinity
	(conf.GetProperty<string>(RtConfig::RtCpuAffinity), affinity);
      size_t nr_threads = conf.GetProperty<size_t>(RtConfig::RtNrThreads);

      LOG_INFO << "Number of threads: " << nr_threads << "\n";
      stringstream os;
      os << "Thread affinity: {";
      for (size_t i = 0; i < nr_threads; ++i) {
        if (i != 0)
	  os << ",";
        if (affinity.size() != nr_threads)
	  os << i;
        else
	  os << affinity[i];
      }

      os << "}\n";
      LOG_INFO << os.str();
    }

    void RtCtx::SetRtCtx(const RtConfig &conf)
    {
      size_t nr_threads = conf.GetProperty<size_t>(RtConfig::RtNrThreads);
      cpu_affinity_ = ParseCpuAffinity
	(conf.GetProperty<string>(RtConfig::RtCpuAffinity), cpu_affinity_);

      if (cpu_affinity_.size() != nr_threads) {
        for (size_t i = 1; i < nr_threads; ++i)
	  cpu_affinity_.push_back(i);
      }
    }

  } // end of namespace runtime
} // end of namespace sparsex
