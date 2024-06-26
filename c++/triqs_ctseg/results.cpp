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
// Authors: Nikita Kavokine, Hao Lu, Olivier Parcollet, Nils Wentzell

#include "./results.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "G_tau", c.G_tau);
    h5_write(grp, "sign", c.sign);
    h5_write(grp, "F_tau", c.F_tau);
    h5_write(grp, "nn_tau", c.nn_tau);
    h5_write(grp, "sperp_tau", c.sperp_tau);
    h5_write(grp, "nn_static", c.nn_static);
    h5_write(grp, "densities", c.densities);
    h5_write(grp, "pert_order_histo_Delta", c.pert_order_histo_Delta);
    h5_write(grp, "pert_order_histo_Jperp", c.pert_order_histo_Jperp);
    h5_write(grp, "state_hist", c.state_hist);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, results_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "G_tau", c.G_tau);
    h5_read(grp, "sign", c.sign);
    h5_read(grp, "F_tau", c.F_tau);
    h5_read(grp, "nn_tau", c.nn_tau);
    h5_read(grp, "sperp_tau", c.sperp_tau);
    h5_read(grp, "nn_static", c.nn_static);
    h5_read(grp, "densities", c.densities);
    h5_read(grp, "pert_order_histo_Delta", c.pert_order_histo_Delta);
    h5_read(grp, "pert_order_histo_Jperp", c.pert_order_histo_Jperp);
    h5_read(grp, "state_hist", c.state_hist);
  }

} // namespace triqs_ctseg
