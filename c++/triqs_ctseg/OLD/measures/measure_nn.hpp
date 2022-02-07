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
#include "../block_matrix.hpp"

namespace triqs_ctseg {

/// measurement of equal-time density-density correlation function $\langle n_1
/// n_2 \rangle$
/**
 *Measures $$\chi_{a,b}^{\sigma\sigma'}\equiv\langle
 *n_{a\sigma}n_{b\sigma'}\rangle$$ by using the length and overlaps of segments.
 */
struct measure_nn {
  const qmc_parameters *params;
  const configuration *const config;
  /// matrix of segment overlaps; densities on the diagonal
  block_matrix<double> *overlaps;
  matrix<double> nn;
  double Z;

  measure_nn(const qmc_parameters *params_, const configuration *config_,
             block_matrix<double> *overlaps_)
      : params(params_), config(config_), overlaps(overlaps_),
        nn(config_->n_colors(), config_->n_colors()) {
    Z = 0;
    nn() = 0.0;
  }

  /// accumulation
  void accumulate(double s) {
    Z += s;
    for (int i = 0; i < config->n_colors(); i++)
      for (int j = 0; j < config->n_colors(); j++)
        nn(i, j) += s * config->trace.overlap_matrix(i, j);
  }

  /// collect
  void collect_results(mpi::communicator const &c) {
    Z = mpi::all_reduce(Z, c);
    nn = mpi::all_reduce(nn, c);
    nn /= Z * params->beta;

    for (int i = 0; i < config->n_colors(); i++) {
      auto bl_ind_i = config->color_to_block_and_inner_index(i);
      for (int j = 0; j < config->n_colors(); j++) {
        auto bl_ind_j = config->color_to_block_and_inner_index(j);
        (*overlaps)[bl_ind_i.first * config->n_blocks() + bl_ind_j.first](
            bl_ind_i.second, bl_ind_j.second) = nn(i, j);
      }
    }
  }
};
} // namespace triqs_ctseg
