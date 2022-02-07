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

/// Measure the histogram of perturbation orders in $\mathcal{J}_\perp$
/**
 *Measure the histogram of the hybridization perturbation order in
 *$\mathcal{J}_\perp$
 */
struct measure_hist_composite {

  // a histogram
  const int Nmax;
  const configuration *const config;
  matrix<double> *H;
  int N, N_tot;

  // constructor
  measure_hist_composite(configuration const *config_, matrix<double> *H_)
      : Nmax(500), config(config_), H(H_), N(0) {
    H->resize(1, Nmax);
    (*H)() = 0;
  }

  void accumulate(double s) {
    (*H)(0, config->boson_lines.size()) += 1;
    N += 1;
  }

  void collect_results(mpi::communicator const &c) {
    auto H_tot = *H;
    H_tot() = 0.0;

    H_tot = mpi::all_reduce(*H, c);
    N_tot = mpi::all_reduce(N, c);

    for (int i = 0; i < Nmax; i++)
      (*H)(0, i) = H_tot(0, i) / N_tot;
  }
};
} // namespace triqs_ctseg
