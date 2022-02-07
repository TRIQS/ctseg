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
/// Remove a segment
/**
 *Remove a segment $c(\tau)c^\dagger(\tau')$
 */
class move_remove_segment {
  const qmc_parameters *params;
  configuration *config;
  triqs::mc_tools::random_generator &RND;
  int color;

public:
  move_remove_segment(const qmc_parameters *params_, configuration *config_,
                      triqs::mc_tools::random_generator &RND_);

  double attempt();
  double accept();
  void reject();
};
} // namespace triqs_ctseg
