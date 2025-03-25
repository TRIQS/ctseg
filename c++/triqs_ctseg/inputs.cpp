// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "inputs.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, inputs_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "Delta", c.Delta);
    h5_write(grp, "Jperpt", c.Jperpt);
    h5_write(grp, "D0t", c.D0t);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, inputs_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "Delta", c.Delta);
    h5_read(grp, "Jperpt", c.Jperpt);
    h5_read(grp, "D0t", c.D0t);
  }

} // namespace triqs_ctseg
