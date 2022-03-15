#pragma once

// spdlog
#ifdef EXT_DEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF
#endif
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

#define LOG(...) SPDLOG_TRACE(__VA_ARGS__)

// checked always, even in production
#define ALWAYS_EXPECTS(Condition, ErrorMessage, ...)                                                                   \
  if (not(Condition)) {                                                                                                \
    SPDLOG_CRITICAL(ErrorMessage, __VA_ARGS__);                                                                        \
    throw std::runtime_error("Assertion Error, cf log");                                                               \
  }

