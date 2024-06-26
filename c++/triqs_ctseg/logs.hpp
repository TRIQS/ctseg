// Copyright (c) 2022-2024 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Nikita Kavokine, Olivier Parcollet, Nils Wentzell

#pragma once

// Set log levels
#ifdef CTSEG_DEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
static constexpr bool ctseg_debug = true;
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
static constexpr bool ctseg_debug = false;
#endif
#ifdef PRINT_LOGS
static constexpr bool print_logs = true;
#else
static constexpr bool print_logs = false;
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
