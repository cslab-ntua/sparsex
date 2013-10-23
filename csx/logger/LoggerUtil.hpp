#ifndef LOGGER_UTIL_HPP
#define LOGGER_UTIL_HPP

#ifdef __cplusplus
extern "C" {
#endif

void DisableLogging();
void DisableError();
void DisableWarning();
void DisableInfo();
void DisableDebug();
void AlwaysUseConsole();
void AlwaysUseFile(const char *log_file);

void log_error(const char *msg);
void log_warning(const char *msg);

#ifdef __cplusplus
}
#endif

#endif // LOGGER_UTIL_HPP

