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
// Authors: Nikita Kavokine, Hao Lu, Nils Wentzell

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
