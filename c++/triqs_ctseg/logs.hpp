#pragma once

// FIXME : ??
// spdlog
#ifdef PRINT_CONFIG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#elif defined(EXT_DEBUG)
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
//static constexpr bool ctseg_debug = true;
//#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE // FIXME: this doesn't work. Set in work_data
//#else
//static constexpr bool ctseg_debug = false;
//#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF
#endif

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

// Log messages for dubugging
#define LOG(...) SPDLOG_DEBUG(__VA_ARGS__)

// Checked always, even in production
#define ALWAYS_EXPECTS(Condition, ErrorMessage, ...)                                                                   \
  if (not(Condition)) {                                                                                                \
    SPDLOG_CRITICAL(ErrorMessage, __VA_ARGS__);                                                                        \
    throw std::runtime_error("Assertion Error, cf log");                                                               \
  }
