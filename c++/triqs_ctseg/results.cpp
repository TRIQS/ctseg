/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "./results.hpp"

void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

  h5_write(grp, "G_tau", c.G_tau);
  h5_write(grp, "F_tau", c.F_tau);
  h5_write(grp, "K_tau", c.K_tau);
  h5_write(grp, "Kprime_tau", c.Kprime_tau);
  h5_write(grp, "nn_tau", c.nn_tau);
  h5_write(grp, "sperp_tau", c.sperp_tau);
  h5_write(grp, "sperp_tau2", c.sperp_tau2);
  h5_write(grp, "nn_static", c.nn_static);
  h5_write(grp, "densities", c.densities);
  h5_write(grp, "perturbation_order_histo_Delta", c.perturbation_order_histo_Delta);
  h5_write(grp, "perturbation_order_histo_Jperp", c.perturbation_order_histo_Jperp);
}

//------------------------------------

void h5_read(h5::group h5group, std::string subgroup_name, results_t &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

  h5_read(grp, "G_tau", c.G_tau);
  h5_read(grp, "F_tau", c.F_tau);
  h5_read(grp, "K_tau", c.K_tau);
  h5_read(grp, "Kprime_tau", c.Kprime_tau);
  h5_read(grp, "nn_tau", c.nn_tau);
  h5_read(grp, "sperp_tau", c.sperp_tau);
  h5_read(grp, "sperp_tau2", c.sperp_tau2);
  h5_read(grp, "nn_static", c.nn_static);
  h5_read(grp, "densities", c.densities);
  h5_read(grp, "perturbation_order_histo_Delta", c.perturbation_order_histo_Delta);
  h5_read(grp, "perturbation_order_histo_Jperp", c.perturbation_order_histo_Jperp);
}
