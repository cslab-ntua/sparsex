/* -*- C++ -*-
 *
 * Logger_deprecated.hpp --  Logging interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_LOGGER_H__
#define LIBCSX_LOGGER_H__

#include "LogPolicy.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <pthread.h>

using namespace std;

/* The available logging policies */
#ifdef LOGFILE
#define policy FileLogPolicy
#else
#define policy StdLogPolicy
#endif

/* The available logging levels */
namespace logging {
enum LogLevelType {
    error,
    warning,
    info,
    debug
};
}

#ifdef LOGERROR
#define LOG_ERROR(x) LOG<logging::error>(__FILE__,__LINE__,__FUNCTION__, x)
#else
#define LOG_ERROR(...)
#endif

#ifdef LOGWARNING
#define LOG_WARNING(x) LOG<logging::warning>(__FILE__,__LINE__,__FUNCTION__, x)
#else
#define LOG_WARNING(...)
#endif

#ifdef LOGINFO
#define LOG_INFO LogInfo
#else
#define LOG_INFO(...)
#endif

#ifdef LOGDEBUG
#define LOG_DEBUG LogDebug
#else
#define LOG_DEBUG(...)
#endif

/* Generates a unique type for every log_level (error, warning, info, debug) */
template<int T>
struct LogLevel
{
    enum { value = T };
};

template<typename LogPolicy>
class Logger
{
public:
    // In C++11 this is thread safe
    static Logger& GetInstance()
    {
        static Logger obj;
        return obj;
    }

    template<logging::LogLevelType Level, typename... Targs>
    void Print(Targs... args);
    template<logging::LogLevelType Level, typename... Targs>
    void PrintEnhanced(const char *file, unsigned long lineno,
                       const char *function, Targs... args);
private:
    LogPolicy *policy_;
    stringstream log_stream_;
    stringstream postfix_;
    pthread_mutex_t write_mutex_;

    Logger();
    Logger(const Logger &); // Do not implement
    Logger &operator=(const Logger &); // Do not implement

    ~Logger()
    {
        delete policy_;
        pthread_mutex_destroy(&write_mutex_);
    }

    void DoPrint();
    template<typename Arg1, typename... Args>
    void DoPrint(Arg1 arg1, Args... args);

    void AddPrefix(LogLevel<logging::error>)
    {
        log_stream_ << "[ERROR]: ";
    }

    void AddPrefix(LogLevel<logging::warning>)
    {
        log_stream_ << "[WARNING]: ";
    }

    void AddPrefix(LogLevel<logging::info>)
    {   
        log_stream_ << "[INFO]: ";
    }

    void AddPrefix(LogLevel<logging::debug>)
    {
        log_stream_ << "[DEBUG]: ";
    }
};

template<typename LogPolicy>
Logger<LogPolicy>::Logger()
    : policy_(0),
      log_stream_(),
      postfix_()
{
    policy_ = new LogPolicy;
    if (!policy_) {
        cerr << "[LOGGER]: failed to create a Logger instance." 
             << endl;
        exit(1);
    }
    try {
        policy_->open("libcsx_logfile");
    } catch (ios_base::failure &error) {
        cerr << error.what() << endl;
    }
    int err = pthread_mutex_init(&write_mutex_, 0);
    if (err != 0) {
        cerr << "[LOGGER]: failed to initialize lock." 
             << endl;
        exit(1);
    }
}

template<typename LogPolicy>
template<logging::LogLevelType Level, typename... Targs>
void Logger<LogPolicy>::Print(Targs... args)
{
    int err = pthread_mutex_lock(&write_mutex_);
    if (err != 0) {
        cerr << "[LOGGER]: failed to lock shared variable." 
             << endl;
        exit(1);
    }

    AddPrefix(LogLevel<Level>());
    DoPrint(args...);

    err = pthread_mutex_unlock(&write_mutex_);
    if (err != 0) {
        cerr << "[LOGGER]: failed to unlock shared variable." 
             << endl;
        exit(1);
    }
}

template<typename LogPolicy>
template<logging::LogLevelType Level, typename... Targs>
void Logger<LogPolicy>::PrintEnhanced(const char *file, unsigned long lineno,
                                       const char *function, Targs... args)
{
    postfix_ << " [" << file << ":" << lineno << ":" << function << "()]";
    Print<Level>(args...);
}

/* Recursive variadic function */
template<typename LogPolicy>
template<typename Arg1, typename... Args>
void Logger<LogPolicy>::DoPrint(Arg1 arg1, Args... args)
{
    log_stream_ << arg1;
    DoPrint(args...);
}

/* This version will be called when there are no arguments left
   to add to the log_stream_ and the actual write to the log file
   will take place */
template<typename LogPolicy>
void Logger<LogPolicy>::DoPrint()
{
    policy_->write(log_stream_.str() + postfix_.str());
    /* Empty the helper stream */
    log_stream_.str("");
    postfix_.str("");
}

template<typename... Targs>
static inline void LogInfo(Targs... args)
{
    Logger<policy> &logger = Logger<policy>::GetInstance();
    logger.Print<logging::info>(args...);
}

template<typename... Targs>
static inline void LogDebug(Targs... args)
{
    Logger<policy> &logger = Logger<policy>::GetInstance();
    logger.Print<logging::debug>(args...);
}

template<logging::LogLevelType Level, typename... Targs>
void LOG(const char *file, unsigned long lineno, const char *function,
         Targs... args)
{
    Logger<policy> &logger = Logger<policy>::GetInstance();
    logger.PrintEnhanced<Level>(file, lineno, function, args...);
}

#endif // LIBCSX_LOGGER_H__
