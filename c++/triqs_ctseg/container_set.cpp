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
#include "./container_set.hpp"

void h5_write(h5::group h5group, std::string subgroup_name, container_set_t const &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

  h5_write(grp, "G_tau", c.G_tau);
  h5_write(grp, "chi_tau", c.chi_tau);
  h5_write(grp, "densities", c.densities);
}

//------------------------------------

void h5_read(h5::group h5group, std::string subgroup_name, container_set_t &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

  h5_read(grp, "G_tau", c.G_tau);
  h5_read(grp, "chi_tau", c.chi_tau);
  h5_read(grp, "densities", c.densities);
}
