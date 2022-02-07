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
/// Insert spin segment
/**
 *Move which inserts a spin segment $s_+(\tau_1) s_-(\tau_2)$.
 */
class move_insert_spin_segment {

  const qmc_parameters *params;
  configuration *config;
  triqs::mc_tools::random_generator &RND;
  segment seg1;
  bool special_case;

public:
  move_insert_spin_segment(const qmc_parameters *params_,
                           configuration *config_,
                           triqs::mc_tools::random_generator &RND_);
  //----------------------------------------------
  /// from a time, color, compute (length to the left neighbour, length to the
  /// right neighbour, iterator to the right neighbour
  // if no operator, return {beta, beta, end(color)}
  std::tuple<qmc_time_t, qmc_time_t, colored_const_iterator>
  max_lens(qmc_time_t const &tau, int color);
  double attempt();
  double accept();
  void reject();
};
} // namespace triqs_ctseg
