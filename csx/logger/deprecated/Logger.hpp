/* -*- C++ -*-
 *
 * Logger.hpp --  Logging interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_LOGGER_HPP__
#define LIBCSX_LOGGER_HPP__

#include "OutputPolicy.hpp"
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>

namespace logging {

/* The available logging levels */
enum Level {
    error,
    warning,
    info,
    debug
};

template<Level lvl>
struct LevelT
{
    static const char *description() {
        return "[ERROR]: "; //FIXME
    }
};

class SeverityHandler
{
public:
    static bool IsSeverityLogged(Level level)
    {
        return severity_map_[level];
    }

    static void Enable(Level level)
    {
        severity_map_[level] = true;
    }

    static void Disable(Level level)
    {
        severity_map_[level] = false;
    }

private:
    static boost::unordered_map<Level, bool> severity_map_;
};

boost::unordered_map<Level, bool> SeverityHandler::SeverityMap =
                      boost::assign::map_list_of 
                      (error, true)
                      (warning, true)
                      (info, true)
                      (debug, true);


// The template parameter Level is needed to disable logging levels at compile-time
// through template specialization.
template<Level level, typename OutputPolicy = StdPolicy>
class Logger
{
public:
    typedef OutputPolicy output_t;

    static output_t& logging()
    {
        return Singleton<output_t>::GetInstance();
    }

private:
    template<typename T = OutputPolicy>
    struct Singleton
    {
        // In C++11 this is thread safe
        static T& GetInstance()
        {
            static T obj;
            return obj;
        }
    };
};

template<Level level>
static inline
typename Logger<level>::output_t& log() {
    return Logger<level>::logging() << LevelT<level>::description();
}

template<typename T>
static inline
NullPolicy& operator<<(NullPolicy &p, T t) {
    return p;
}

} // end of namespace logging

using namespace logging;

/* Logging macros */
#define LOG_ERROR if(SeverityHandler::IsSeverityLogged(error))  \
        log<error>() << "(" << __FILE__ << ":" << __LINE__      \
                     << ":" << __FUNCTION__ << ") "
#define LOG_WARNING if(SeverityHandler::IsSeverityLogged(warning))  \
        log<warning>() << "(" << __FILE__ << ":" << __LINE__        \
                       << ":" << __FUNCTION__ << ") "
#define LOG_INFO if(SeverityHandler::IsSeverityLogged(info)) log<info>()
#define LOG_DEBUG if(SeverityHandler::IsSeverityLogged(debug)) log<debug>()

/* Compile-time disabling of logging levels */
#define DISABLE_LOGGING_LEVEL(lvl)                                    \
    namespace logging {                                                 \
    template<Level level, typename OutputPolicy> class Logger;       \
    template<typename OutputPolicy> class Logger<lvl, OutputPolicy>;  \
                                                                        \
    template<typename OutputPolicy>                                     \
    class Logger<lvl, OutputPolicy> {                                 \
    public:                                                             \
    typedef NullPolicy output_t;                                        \
    static output_t& logging() {                                        \
        return *reinterpret_cast<output_t *>(0x0);                      \
    }                                                                   \
    };                                                                  \
}

#endif // LIBCSX_LOGGER_HPP__

// #define LOG_ERROR logging::log<logging::Error>()                        \
//     << "(" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ") "
// #define LOG_WARNING logging::log<logging::Warning>()                    \
//     << "(" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ") "
//#define LOG_INFO logging::log<logging::Info>()
//#define LOG_DEBUG logging::log<logging::Debug>()

// typedef boost::function<StdPolicy& ()> LogHandler;

// static LogHandler ErrorHandler = logging::log<logging::Error>;
// static LogHandler WarningHandler = logging::log<logging::Warning>;
// static LogHandler InfoHandler = logging::log<logging::Info>;
// static LogHandler DebugHandler = logging::log<logging::Debug>;

// struct Error
// {
//     static Level level()
//     {
//         return Level::error;
//     }

//     static const char *description() {
//         return "[ERROR]: ";
//     }   
// };

// struct Warning
// {
//     static LevelType level()
//     {
//         return LevelType::warning;
//     }

//     static const char *description() {
//         return "[WARNING]: ";
//     }   
// };

// struct Info
// {
//     static LevelType level()
//     {
//         return LevelType::info;
//     }

//     static const char *description() {
//         return "[INFO]: ";
//     }   
// };

// struct Debug
// {
//     static LevelType level()
//     {
//         return LevelType::debug;
//     }

//     static const char *description() {
//         return "[DEBUG]: ";
//     }
// };
