/*
 * LoggerUtil.hpp --  Logging utilities.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LOGGER_UTIL_HPP
#define LOGGER_UTIL_HPP

#include "../cdecl.h"

BEGIN_C_DECLS

void DisableLogging();
void DisableError();
void DisableWarning();
void DisableInfo();
void DisableDebug();
void AlwaysUseConsole();
void AlwaysUseFile(const char *log_file);

void log_error(const char *msg);
void log_warning(const char *msg);

END_C_DECLS

#endif // LOGGER_UTIL_HPP

