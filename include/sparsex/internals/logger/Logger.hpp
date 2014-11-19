/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Logger.hpp
 * \brief Logging interface
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_LOGGER_LOGGER_HPP
#define SPARSEX_INTERNALS_LOGGER_LOGGER_HPP

#include <sparsex/internals/logger/OutputPolicy.hpp>
#include <sstream>
#include <cstdlib>
#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>

namespace sparsex {
namespace logging {

/* The logging levels */
enum Level {
    Error = 0,
    Warning = 1,
    Info = 2,
    Verbose = 3,
    Debug = 4,
    None = 5
};

/* The available sinks */
enum Sink {
    Null,
    Console,
    File
};

#define DEFAULT_ERROR_POLICY Console
#define DEFAULT_WARNING_POLICY Console
#define DEFAULT_INFO_POLICY Null
#define DEFAULT_VERBOSE_POLICY Null
#define DEFAULT_DEBUG_POLICY Null

class LoggingHandler
{
    typedef boost::unordered_map<Level, OutputHandler>::iterator iterator;

public:
    static LoggingHandler &GetInstance()
    {
        static LoggingHandler instance;
        return instance;
    }

    OutputHandler GetHandler(Level level)
    {
        return handlers_[level];
    }

    void SetHandler(Level level, Sink sink)
    {
        handlers_[level] = boost::bind(boost::ref(sinks_[sink]), _1);
    }

    void SetHandlers(Level level, Sink sink);
    void SetLogFile(const char *log_file);

    static string &Prefix(Level level)
    {
        return names_[level];
    }

private:
    static boost::unordered_map<Level, string> names_;
    boost::unordered_map<Sink, OutputHandler> sinks_;
    boost::unordered_map<Level, OutputHandler> handlers_;

    LoggingHandler();
    ~LoggingHandler() {}
    LoggingHandler(const LoggingHandler &); // Do not implement
    LoggingHandler &operator=(const LoggingHandler &); // Do not implement
};

/* Proper logging with this class is based on the fact that:
   - Temporary objects are destroyed as the last step in evaluating the 
     full-expression that (lexically) contains the point where they were 
     created.
   This allows for the actual writing to the output sink to be performed 
   in the dtor. */
template<Level level>
class Logger
{
public:
    ~Logger()
    {
        LoggingHandler& h = LoggingHandler::GetInstance();
        h.GetHandler(level)(stream_buffer_);
    }

    template<typename T>
    Logger& operator<<(const T &t)
    {
        stream_buffer_ << t;
        return *this;
    }

private:
    ostringstream stream_buffer_;
};

template<>
class Logger<None>
{
public:
    ~Logger() {}

    template<typename T>
    Logger& operator<<(const T &t)
    {
        return *this;
    }
};

void DisableLevel(Level level);
void EnableConsole(Level level);
void EnableFile(Level level, const char *log_file);

} // end of namespace logging
} // end of namespace sparsex

using namespace sparsex::logging;

/* Logging macros */
#define LOG_ERROR Logger<Error>() << LoggingHandler::Prefix(Error) \
        << "(" << __FILE__                                          \
        << ":" << __LINE__                                          \
        << ":" << __FUNCTION__ << "()) "

#define LOG_WARNING Logger<Warning>() << LoggingHandler::Prefix(Warning) \
        << "(" << __FILE__                                              \
        << ":" << __LINE__                                              \
        << ":" << __FUNCTION__ << "()) "

#define LOG_INFO Logger<Info>() << LoggingHandler::Prefix(Info)
#define LOG_VERBOSE Logger<Verbose>() << LoggingHandler::Prefix(Verbose)

#if SPX_DEBUG
#   define LOG_DEBUG Logger<Debug>() << LoggingHandler::Prefix(Debug)
#else
#   define LOG_DEBUG Logger<None>()
#endif

#endif // SPARSEX_INTERNALS_LOGGER_LOGGER_HPP
