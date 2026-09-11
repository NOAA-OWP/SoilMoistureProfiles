#ifndef SMP_LOGGER_HPP
#define SMP_LOGGER_HPP

#ifdef SMP_USE_EWTS
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Using EWTS
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#include "ewts/module_constants.hpp"
#include "ewts/logger.hpp"
#include "ewts/log_levels.hpp"

#define LOG(...) ::ewts::GetLogger(::ewts::modules::EWTS_ID_SMP).Log(__VA_ARGS__)
#define GetLogLevel() ::ewts::GetLogger(::ewts::modules::EWTS_ID_SMP).GetLogLevel()
#define IsLoggingEnabled() ::ewts::GetLogger(::ewts::modules::EWTS_ID_SMP).IsLoggingEnabled()

using ewts::EwtsInit;
using ewts::LogLevel;

inline constexpr const char* SMP_MODULE_ID = ewts::modules::EWTS_ID_SMP;

#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Log messages written to STDOUT
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#include <chrono>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>

inline constexpr const char* SMP_MODULE_ID = "SMP";

// Minimal LogLevel enum (match EWTS values if possible)
enum class LogLevel {
    NOTSET = 0,
    DEBUG = 10,
    INFO = 20,
    WARNING = 30,
    SEVERE = 40,
    FATAL = 50
};

#define SMP_FALLBACK_LOGGING_ENABLED   1
#define SMP_FALLBACK_LOG_LEVEL         LogLevel::INFO

inline std::string smp_utc_timestamp() {
    using clock = std::chrono::system_clock;
    auto now = clock::now();
    auto tt = clock::to_time_t(now);

    std::tm tm{};
    gmtime_r(&tt, &tm);

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S")
        << '.' << std::setw(3) << std::setfill('0') << ms.count()
        << 'Z';

    return oss.str();
}

inline std::string_view smp_rstrip_newline(std::string_view s) {
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) {
        s.remove_suffix(1);
    }
    return s;
}

inline const char* smp_level_to_string(LogLevel level) {
    switch (level) {
        case LogLevel::DEBUG:   return "DEBUG";
        case LogLevel::INFO:    return "INFO";
        case LogLevel::WARNING: return "WARN";
        case LogLevel::SEVERE:  return "ERROR";
        case LogLevel::FATAL:   return "FATAL";
        default:                return "LOG";
    }
}

inline bool IsLoggingEnabled() {
    return SMP_FALLBACK_LOGGING_ENABLED;
}

inline LogLevel GetLogLevel() {
    return SMP_FALLBACK_LOG_LEVEL;  // simple default
}

inline void Log(LogLevel level, std::string_view message) {

    if (message.empty()) return;
    if (!IsLoggingEnabled()) return;
    if (level < GetLogLevel()) return;

    std::ostringstream oss;

    message = smp_rstrip_newline(message);

    oss << smp_utc_timestamp() << " "
        << SMP_MODULE_ID << " "
        << smp_level_to_string(level) << " "
        << message << '\n';

    std::cout << oss.str() << std::flush;   
}

inline void Log(std::string_view message, LogLevel level) {
    Log(level, message);
}

inline void Log(LogLevel level, const char* fmt, ...) {
    if (!fmt) return;
    if (!IsLoggingEnabled()) return;
    if (level < GetLogLevel()) return;

    va_list args;
    va_start(args, fmt);

    va_list args_copy;
    va_copy(args_copy, args);
    int len = std::vsnprintf(nullptr, 0, fmt, args_copy);
    va_end(args_copy);

    if (len < 0) {
        va_end(args);
        return;
    }

    std::string buffer(static_cast<std::size_t>(len) + 1, '\0');
    std::vsnprintf(buffer.data(), buffer.size(), fmt, args);
    va_end(args);

    buffer.resize(static_cast<std::size_t>(len));

    Log(level, buffer);
}

#define LOG(...) Log(__VA_ARGS__)
#endif
#endif /* SMP_LOGGER_HPP */
