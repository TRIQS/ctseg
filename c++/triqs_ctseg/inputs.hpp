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

#include "params.hpp"
#include <triqs/gfs.hpp>

namespace triqs_ctseg {

  // Group the inputs of the solver.

  struct inputs_t {

    block_gf<imtime> Delta;           // hybridization function
    gf<imtime, matrix_valued> Jperpt; // perpendicular spin-spin interaction
    block2_gf<imtime> D0t;            // retarded density-density interaction
  };

  // h5_read/write
  void h5_write(h5::group h5group, std::string subgroup_name, inputs_t const &s);
  void h5_read(h5::group h5group, std::string subgroup_name, inputs_t &s);

} // namespace triqs_ctseg
