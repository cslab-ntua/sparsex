/* -*- C++ -*-
 *
 * logger.cc --  Logging interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "logger.hpp"

namespace logging {

LoggingHandler::LoggingHandler()
{
    sinks_[none] = boost::ref(NullPolicy::GetInstance());
    sinks_[console] = boost::ref(StdPolicy::GetInstance());
    sinks_[file] = boost::ref(FilePolicy::GetInstance());

    /* Set default policy for each level */
    handlers_[error] = boost::ref(sinks_[DEFAULT_ERROR_POLICY]);
    handlers_[warning] = boost::ref(sinks_[DEFAULT_WARNING_POLICY]);
    handlers_[info] = boost::ref(sinks_[DEFAULT_INFO_POLICY]);
    handlers_[debug] = boost::ref(sinks_[DEFAULT_DEBUG_POLICY]);
}

void LoggingHandler::SetHandlers(Sink sink)
{
    for (iterator it = handlers_.begin(); it != handlers_.end(); ++it)
        (*it).second = boost::bind(boost::ref(sinks_[sink]), _1);
}

void DisableLogging()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(none);
}

void DisableLevel(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(level, none);
}

void UseConsole(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(level, console);
}

void AlwaysUseConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(console);
}

void UseFile(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile("logfile");
    h.SetHandler(level, file);
}

void UseFile(Level level, const char *log_file)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile(log_file);
    h.SetHandler(level, file);
}

void AlwaysUseFile()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile("logfile");
    h.SetHandlers(file);
}

void AlwaysUseFile(const char *log_file)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile(log_file);
    h.SetHandlers(file);
}

}
