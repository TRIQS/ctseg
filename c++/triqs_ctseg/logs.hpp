#pragma once

// Set log levels
#ifdef CTSEG_DEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
static constexpr bool ctseg_debug = true;
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
static constexpr bool ctseg_debug      = false;
#endif
#ifdef CHECK_INVARIANTS
static constexpr bool check_invariants = true;
#else
static constexpr bool check_invariants = false;
#endif

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

// Log messages for dubugging
#define LOG(...) SPDLOG_DEBUG(__VA_ARGS__)

// Checked always, even in production
#define ALWAYS_EXPECTS(Condition,  ...)                                                                   \
  if (not(Condition)) {                                                                                                \
    SPDLOG_CRITICAL("Error in function {} in file  {} at line {}", __FUNCTION__, __FILE__, __LINE__);                  \
    SPDLOG_CRITICAL( __VA_ARGS__);                                                                      \
    throw std::runtime_error("Assertion Error, cf log");                                                               \
  }
