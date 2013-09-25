/* -*- C++ -*-
 *
 * logger.hpp --  Logging interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_LOGGER_HPP__
#define LIBCSX_LOGGER_HPP__

#include "output_policy.hpp"

#include <sstream>
#include <cstdlib>
#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>

namespace logging {

/* The logging levels */
enum Level {
    error,
    warning,
    info,
    debug
};

/* The available sinks */
enum Sink {
    none,
    console,
    file
};

#define DEFAULT_ERROR_POLICY console
#define DEFAULT_WARNING_POLICY console
#define DEFAULT_INFO_POLICY console
#define DEFAULT_DEBUG_POLICY none

void DisableLogging();
void DisableLevel(Level level);
void UseConsole(Level level);
void AlwaysUseConsole();
void UseFile(Level level);
void UseFile(Level level, const char *log_file);
void AlwaysUseFile();
void AlwaysUseFile(const char *log_file);

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

    void SetHandlers(Sink sink);

    void SetLogFile(const char *log_file)
    {
        FilePolicy& f = FilePolicy::GetInstance();
        f.Init(log_file);
    }

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

boost::unordered_map<Level, string> LoggingHandler::names_ =
    boost::assign::map_list_of 
    (error, "[ERROR]: ")
    (warning, "[WARNING]: ")
    (info, "[INFO]: ")
    (debug, "[DEBUG]: ");

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

} // end of namespace logging

using namespace logging;

/* Logging macros */
#define LOG_ERROR Logger<error>() << LoggingHandler::Prefix(error) \
        << "(" << __FILE__                                          \
        << ":" << __LINE__                                          \
        << ":" << __FUNCTION__ << "()) "
#define LOG_WARNING Logger<warning>() << LoggingHandler::Prefix(warning) \
        << "(" << __FILE__                                              \
        << ":" << __LINE__                                              \
        << ":" << __FUNCTION__ << "()) "
#define LOG_INFO Logger<info>() << LoggingHandler::Prefix(info)
#define LOG_DEBUG Logger<debug>() << LoggingHandler::Prefix(debug)

#endif // LIBCSX_LOGGER_HPP__
