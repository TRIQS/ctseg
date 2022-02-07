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

constexpr unsigned int ipow(unsigned int n, unsigned int m) { // move to
                                                              // utility?
  return m == 0 ? 1 : m == 1 ? n : n * ipow(n, m - 1);
}

namespace triqs_ctseg {

/// Measure the state histogram
/**
 *Measures the time the impurity spends in a certain atomic state:
 * an atomic state is characterized by occupation number for each color, e.g.
 * $|n_1 n_2 n_3 n_4\rangle = |0 1 1 0\rangle$ for a model with 4 colors
 *
 * - the index of a state in the histogram is given by $\sum_i n_i 2^i$
 *
 * - the length of the histogram is 2^n_colors
 */

struct measure_statehist {

  const configuration *config;
  const qmc_parameters *params;
  vector<double> *H;
  int N, N_tot;
  const double beta;

  // constructor
  measure_statehist(const qmc_parameters *params_, const configuration *config_,
                    vector<double> *H_)
      : params(params_), config(config_), H(H_), N(0), beta(params->beta) {
    H->resize(ipow(2, config->n_colors()));
    (*H)() = 0;
  }

  void accumulate(double s) {

    N += 1;

    double tau_prev = 0.0; // time of prevous operator; start with 0
    for (auto it = config->ops_map().begin(); it != config->ops_map().end();
         ++it) {
      int state = 0; // get state for this operator
      for (size_t col = 0; col < config->n_colors(); ++col) {
        // get the index of the impurity state
        if (it->right_occupation()[col])
          state += ipow(2, col);
      }
      (*H)(state) += (it->tau - tau_prev); // increment
      tau_prev = (double)it->tau;
    }

    int state = 0;
    // get edge state contribution; tau_prev has time of last operator
    for (size_t col = 0; col < config->n_colors(); ++col) {
      if (config->edge_state(col))
        state += ipow(2, col);
    }
    (*H)(state) += (beta - tau_prev); // increment
  }

  void collect_results(mpi::communicator const &c) {
    auto H_tot = *H;
    H_tot() = 0.0;

    H_tot = mpi::all_reduce(*H, c);
    N_tot = mpi::all_reduce(N, c);

    (*H) = H_tot / (N_tot * beta);
  }
};
} // namespace triqs_ctseg
