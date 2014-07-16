/*
 * Copyright (C) 2013-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Logger.cpp
 * \brief Logging interface
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/logger/Logger.hpp>
#include <sparsex/internals/logger/LoggerUtil.hpp>

namespace sparsex {
namespace logging {

boost::unordered_map<Level, string> LoggingHandler::names_ =
    boost::assign::map_list_of 
    (Error, "[ERROR]: ")
    (Warning, "[WARNING]: ")
    (Info, "[INFO]: ")
    (Verbose, "[VERBOSE]: ")
    (Debug, "[DEBUG]: ");

LoggingHandler::LoggingHandler()
{
    sinks_[Null] = boost::ref(NullPolicy::GetInstance());
    sinks_[Console] = boost::ref(StdPolicy::GetInstance());
    sinks_[File] = boost::ref(FilePolicy::GetInstance());

    /* Set default policy for each level */
    handlers_[Error] = boost::ref(sinks_[DEFAULT_ERROR_POLICY]);
    handlers_[Warning] = boost::ref(sinks_[DEFAULT_WARNING_POLICY]);
    handlers_[Info] = boost::ref(sinks_[DEFAULT_INFO_POLICY]);
    handlers_[Verbose] = boost::ref(sinks_[DEFAULT_VERBOSE_POLICY]);
    handlers_[Debug] = boost::ref(sinks_[DEFAULT_DEBUG_POLICY]);
}

void LoggingHandler::SetHandlers(Level level, Sink sink)
{
    for (iterator it = handlers_.begin(); it != handlers_.end(); ++it)
        if ((*it).first <= level)
            (*it).second = boost::bind(boost::ref(sinks_[sink]), _1);
}

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

void EnableConsole(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(level, Console);
}

void EnableFile(Level level)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile("logfile");
    h.SetHandlers(level, File);
}

} // end of namespace logging
} // end of namespace sparsex

void UseFile(const char *log_file)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile(log_file);
}

void DisableLogging()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(None, Null);
}

void EnableAllConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Debug, Console);
}

void EnableAllFile(const char *log_file)
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetLogFile(log_file);
    h.SetHandlers(Debug, File);
}

void EnableErrorConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Error, Console);
}

void EnableErrorFile()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Error, File);
}

void EnableWarningConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Warning, Console);
}

void EnableWarningFile()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Warning, File);
}

void EnableInfoConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Info, Console);
}

void EnableInfoFile()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Info, File);
}

void EnableVerboseConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Verbose, Console);
}

void EnableVerboseFile()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Verbose, File);
}

void EnableDebugConsole()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Debug, Console);
}

void EnableDebugFile()
{
    LoggingHandler& h = LoggingHandler::GetInstance();
    h.SetHandlers(Debug, File);
}

void log_error(const char *msg)
{
    Logger<Error>() << LoggingHandler::Prefix(Error) << msg;
}

void log_warning(const char *msg)
{
    Logger<Warning>() << LoggingHandler::Prefix(Warning) << msg;
}

void log_verbose(const char *msg)
{
    Logger<Verbose>() << LoggingHandler::Prefix(Verbose) << msg;
}
