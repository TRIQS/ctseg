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

namespace triqs_ctseg {

/// measurement of the MC sign
/**
  *Measures the Monte-Carlo sign
  @note This sign is 1 except in the presence of off-diagonal hybridization or
  positive $\mathcal{J}_\perp(\tau)$.
  */
struct measure_sign {
  double *Z;
  int N;

  measure_sign(double *Z_) : Z(Z_), N(0) { (*Z) = 0.0; }

  /// accumulation
  void accumulate(double s) {
    (*Z) += s;
    N++;
  }

  /// collect
  void collect_results(mpi::communicator const &c) {
    N = mpi::all_reduce(N, c);
    (*Z) = mpi::all_reduce((*Z), c);
    (*Z) = (*Z) / N;
  }
};
} // namespace triqs_ctseg
