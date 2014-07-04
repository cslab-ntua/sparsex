/*
 * LoggerUtil.hpp --  Logging utilities.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_LOGGER_LOGGER_UTIL_HPP
#define SPARSEX_INTERNALS_LOGGER_LOGGER_UTIL_HPP

#include "../cdecl.h"

SPX_BEGIN_C_DECLS__

void DisableLogging();
void DisableError();
void DisableWarning();
void DisableInfo();
void DisableDebug();
void AlwaysUseConsole();
void AlwaysUseFile(const char *log_file);

void log_error(const char *msg);
void log_warning(const char *msg);

SPX_END_C_DECLS__

#endif  // SPARSEX_INTERNALS_LOGGER_LOGGER_UTIL_HPP

