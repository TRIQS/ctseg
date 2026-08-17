// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./results.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5::write_hdf5_format(grp, c);
    h5_write(grp, "G_tau", c.G_tau);
    h5_write(grp, "average_sign", c.average_sign);
    h5_write(grp, "F_tau", c.F_tau);
    h5_write(grp, "nn_tau", c.nn_tau);
    h5_write(grp, "nn_nu_dlr", c.nn_nu_dlr);
    h5_write(grp, "Sperp_tau", c.Sperp_tau);
    h5_write(grp, "Sminus_Splus_tau", c.Sminus_Splus_tau);
    h5_write(grp, "Splus_Sminus_tau", c.Splus_Sminus_tau);
    h5_write(grp, "nn_static", c.nn_static);
    h5_write(grp, "densities", c.densities);
    h5_write(grp, "pert_order_Delta", c.pert_order_Delta);
    h5_write(grp, "average_order_Delta", c.average_order_Delta);
    h5_write(grp, "pert_order_Jperp", c.pert_order_Jperp);
    h5_write(grp, "average_order_Jperp", c.average_order_Jperp);
    h5_write(grp, "state_hist", c.state_hist);
    h5_write(grp, "dyn_phi_n", c.dyn_phi_n);
    h5_write(grp, "dyn_phi_phi", c.dyn_phi_phi);
    h5_write(grp, "g2w", c.g2w);
    h5_write(grp, "g3w", c.g3w);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, results_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "G_tau", c.G_tau);
    h5_read(grp, "average_sign", c.average_sign);
    h5_read(grp, "F_tau", c.F_tau);
    h5_read(grp, "nn_tau", c.nn_tau);
    h5_read(grp, "nn_nu_dlr", c.nn_nu_dlr);
    h5_read(grp, "Sperp_tau", c.Sperp_tau);
    h5::try_read(grp, "Sminus_Splus_tau", c.Sminus_Splus_tau);
    h5::try_read(grp, "Splus_Sminus_tau", c.Splus_Sminus_tau);
    h5_read(grp, "nn_static", c.nn_static);
    h5_read(grp, "densities", c.densities);
    h5_read(grp, "pert_order_Delta", c.pert_order_Delta);
    h5_read(grp, "average_order_Delta", c.average_order_Delta);
    h5_read(grp, "pert_order_Jperp", c.pert_order_Jperp);
    h5_read(grp, "average_order_Jperp", c.average_order_Jperp);
    h5_read(grp, "state_hist", c.state_hist);
    h5::try_read(grp, "dyn_phi_n", c.dyn_phi_n);
    h5::try_read(grp, "dyn_phi_phi", c.dyn_phi_phi);
    h5_read(grp, "g2w", c.g2w);
    h5_read(grp, "g3w", c.g3w);
  }

} // namespace triqs_ctseg
