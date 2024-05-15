#pragma once

// Set log levels
#ifdef CTSEG-J_DEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
static constexpr bool ctseg_debug = true;
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
static constexpr bool ctseg_debug = false;
#endif
#ifdef PRINT_LOGS
static constexpr bool print_logs = true;
#else
static constexpr bool print_logs  = false;
#endif

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>

#include <nda/nda.hpp>
template <nda::MemoryArray A> struct fmt::formatter<A> : ostream_formatter {};

// Log messages for dubugging
#define LOG(...) SPDLOG_DEBUG(__VA_ARGS__)

// Checked always, even in production
#define ALWAYS_EXPECTS(Condition, ...)                                                                                 \
  if (not(Condition)) {                                                                                                \
    SPDLOG_CRITICAL("Error in function {} in file  {} at line {}", __FUNCTION__, __FILE__, __LINE__);                  \
    SPDLOG_CRITICAL(__VA_ARGS__);                                                                                      \
    throw std::runtime_error("Assertion Error, cf log");                                                               \
  }
