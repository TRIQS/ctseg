/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once
#include "../configuration.hpp"
#include "../qmc_parameters.hpp"
namespace triqs_ctseg {
/**
  * Move that takes a spin segment and splits it in two separate segments
  @warning Experimental feature
*/
class move_split_spin_segment2 {
  qmc_parameters *params;
  configuration *config;
  triqs::mc_tools::random_generator &RND;
  colored_const_iterator A_up, B_dn, A_dn, B_up, A_dn_new, B_up_new;
  segment seg_not_shifted;

  qmc_time_t B_up_tau, A_dn_tau;
  double length;
  bool tried_split;

public:
  move_split_spin_segment2(qmc_parameters *params_, configuration *config_,
                           triqs::mc_tools::random_generator &RND_);
  double attempt();
  double accept();
  void reject();
};
} // namespace triqs_ctseg
