// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

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
