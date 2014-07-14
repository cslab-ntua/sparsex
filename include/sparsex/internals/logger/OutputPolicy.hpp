/*
 * \file OutputPolicy.hpp
 *
 * \brief Output policies for logging
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_LOGGER_OUTPUT_POLICY_HPP
#define SPARSEX_INTERNALS_LOGGER_OUTPUT_POLICY_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <pthread.h>
#include <boost/function.hpp>

using namespace std;

namespace sparsex {
namespace logging {

template<typename T>
class Singleton
{
public:
    // In C++11 this is thread safe
    static T& GetInstance()
    {
        static T obj;
        return obj;
    }
};

class NullPolicy : public Singleton<NullPolicy>
{
public:
    typedef ostringstream buffer_type;
    void operator()(const buffer_type &buffer) {}
};

/* Logging policy that directs all the output of the logger to stdout */ 
class StdPolicy : public Singleton<StdPolicy>
{
public:
    typedef ostringstream buffer_type;
    void operator()(const buffer_type &buffer);

private:
    pthread_mutex_t write_mutex_;
};

/* Logging policy that directs all the output of the logger to a file */ 
class FilePolicy : public Singleton<FilePolicy>
{
public:
    typedef ostringstream buffer_type;
    void operator()(const buffer_type &buffer);
    void Init(const char *log_file);
    ~FilePolicy()
    {
        os_.close();
    }

private:
    ofstream os_;
    pthread_mutex_t write_mutex_;
};

typedef boost::function<void (ostringstream&)> OutputHandler;

}   // end of namespace logging
}   // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_LOGGER_OUTPUT_POLICY_HPP
