// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <vector>
#include "configuration.hpp"
#include "work_data.hpp"

namespace triqs_ctseg {

  void check_invariant(configuration_t const &config, work_data_t const &wdata);

  void check_segments(configuration_t const &config);

  void check_dets(configuration_t const &config, work_data_t const &wdata);

  void check_jlines(configuration_t const &config);

} // namespace triqs_ctseg
