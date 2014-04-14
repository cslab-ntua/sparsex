/*
 * Logger.cpp --  Logging interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "Logger.hpp"
#include "LoggerUtil.hpp"

namespace logging {

LoggingHandler::LoggingHandler()
{
    sinks_[Null] = boost::ref(NullPolicy::GetInstance());
    sinks_[Console] = boost::ref(StdPolicy::GetInstance());
    sinks_[File] = boost::ref(FilePolicy::GetInstance());

    /* Set default policy for each level */
    handlers_[Error] = boost::ref(sinks_[DEFAULT_ERROR_POLICY]);
    handlers_[Warning] = boost::ref(sinks_[DEFAULT_WARNING_POLICY]);
    handlers_[Info] = boost::ref(sinks_[DEFAULT_INFO_POLICY]);
    handlers_[Debug] = boost::ref(sinks_[DEFAULT_DEBUG_POLICY]);
}

void LoggingHandler::SetHandlers(Sink sink)
{
    for (iterator it = handlers_.begin(); it != handlers_.end(); ++it)
        (*it).second = boost::bind(boost::ref(sinks_[sink]), _1);
}

} // end of namespace logging

void LoggingHandler::SetLogFile(const char *log_file)
{
    FilePolicy& f = FilePolicy::GetInstance();
    f.Init(log_file);
}

void DisableLevel(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(level, Null);
}

void UseConsole(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(level, Console);
}

void UseFile(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile("logfile");
    h.SetHandler(level, File);
}

void UseFile(Level level, const char *log_file)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile(log_file);
    h.SetHandler(level, File);
}

void DisableLogging()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Null);
}

void DisableError()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(Error, Null);
}

void DisableWarning()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(Warning, Null);
}

void DisableInfo()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(Info, Null);
}

void DisableDebug()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandler(Debug, Null);
}

void AlwaysUseConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Console);
}

void AlwaysUseFile(const char *log_file)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile(log_file);
    h.SetHandlers(File);
}

void log_error(const char *msg)
{
    Logger<Error>() << LoggingHandler::Prefix(Error) << msg;
}

void log_warning(const char *msg)
{
    Logger<Warning>() << LoggingHandler::Prefix(Warning) << msg;
}
