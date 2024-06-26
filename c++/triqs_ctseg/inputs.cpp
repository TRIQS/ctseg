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
// Authors: Olivier Parcollet, Nils Wentzell

#include "inputs.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, inputs_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "delta", c.delta);
    h5_write(grp, "jperpt", c.jperpt);
    h5_write(grp, "d0t", c.d0t);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, inputs_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "delta", c.delta);
    h5_read(grp, "jperpt", c.jperpt);
    h5_read(grp, "d0t", c.d0t);
  }

} // namespace triqs_ctseg
