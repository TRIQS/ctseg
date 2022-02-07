/*******************************************************************************
 *
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
 *
 ******************************************************************************/
#pragma once
#include "../configuration.hpp"
#include "../precompute_nw.hpp"

namespace triqs_ctseg {

// measurement for the density-density correlation functionon Matsubara
// frequencies measures <n_1n_2>(omega) as Fourier transform of <n_1(tau)n_2(0)>
// and hence does not exploit time translational invariance

struct measure_nnw {

  // pointer to the MC configuration
  const qmc_parameters *params;
  const configuration *config;

  gf<imfreq> &nnw; //~ Norb*Norb matrix
  const double beta;
  const int n_imfreq; // number of bosonic Matsubara frequencies
  double Z;

  precompute_nw nw;

  // constructor
  measure_nnw(const qmc_parameters *params_, const configuration *config_,
              gf<imfreq> &nnw_)
      : params(params_), config(config_), nnw(nnw_), beta(params->beta),
        n_imfreq(nnw.mesh().size()), Z(0.0), nw(params, config) {
    nnw.data()() = 0;
    nw.minimum_size(n_imfreq);
  }

  void accumulate(double s) {

    Z += s;

    nw(); // precompute nw

    for (int orb2 = 0; orb2 < config->n_colors(); orb2++)
      if (config->edge_state(orb2) == 1) // if n_2(tau=0)=1
        for (int orb1 = 0; orb1 < config->n_colors(); orb1++)
          for (int m = 0; m < n_imfreq; ++m)
            nnw.data()(m, orb1, orb2) += nw(orb1, m);

  } // accumulate

  // average results and normalize nnw
  void collect_results(mpi::communicator const &c) {
    Z   = mpi::reduce(Z, c);
    nnw = mpi::reduce(nnw, c);
    nnw = nnw / Z;
  }
};

} // namespace triqs_ctseg
