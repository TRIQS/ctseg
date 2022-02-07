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
/// Group two segments into a spin segment
/**
  *Move that takes two separate segments and groups them in a spin segment
  *
  @warning Experimental feature
 */
class move_group_into_spin_segment {
  qmc_parameters *params;
  configuration *config;
  triqs::mc_tools::random_generator &RND;

  colored_const_iterator O1, O2, G1p, G2p;
  qmc_time_t O1_tau, O2_tau;
  bool tried_group;

public:
  move_group_into_spin_segment(qmc_parameters *params_, configuration *config_,
                               triqs::mc_tools::random_generator &RND_);
  std::pair<bool, colored_const_iterator>
  find_point_to_group_with(colored_const_iterator O, colored_const_iterator L,
                           colored_const_iterator R);
  double attempt();
  double accept();
  void reject();
};
} // namespace triqs_ctseg
