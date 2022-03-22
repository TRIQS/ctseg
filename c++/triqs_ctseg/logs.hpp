#pragma once

// spdlog
#ifdef EXT_DEBUG
//#include "spdlog/common.h"
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE // FIXME: this doesn't work. Set in work_data
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF
#endif
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

#define LOG(...) SPDLOG_DEBUG(__VA_ARGS__)

// checked always, even in production
#define ALWAYS_EXPECTS(Condition, ErrorMessage, ...)                                                                   \
  if (not(Condition)) {                                                                                                \
    SPDLOG_CRITICAL(ErrorMessage, __VA_ARGS__);                                                                        \
    throw std::runtime_error("Assertion Error, cf log");                                                               \
  }
